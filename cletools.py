import pysam
import pandas as pd
import re, sys

def revcomp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}  # DNA complement pairs
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna))
    return reverse_complement

def convert_aa(codon):
    three = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"]
    one  = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    for i in range(0,len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i],codon)
    return codon

class SmallVariantVcf:
    def __init__(self, file=None,transcripts=None):

        # initialize vcf object
        self.vcf = pysam.VariantFile(file,"r") if isinstance(file, str) else None
        self.__process_transcript_file__(transcripts)

        # initialize dataframe with small variants
        self.variants = pd.DataFrame(columns=['type','filters','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf'])

        self.df = pd.DataFrame()
        self.vep_fields = None

        if self.vcf != None:
            self.parse()

    def __process_transcript_file__(self,transcripts):
        if isinstance(transcripts, pd.DataFrame) and set(['Gene','Transcript','cdsStart','cdsEnd','Strand']).issubset(transcripts.columns):
            transcripts = transcripts.drop(columns='Gene')
            self.transcripts = pd.concat([transcripts[(transcripts['Strand']==1)].groupby(['Transcript','Strand']).agg({'cdsStart': 'min', 'cdsEnd': 'max'}).reset_index(),transcripts[(transcripts['Strand']==-1)].groupby(['Transcript','Strand']).agg({'cdsStart': 'max', 'cdsEnd': 'min'}).reset_index()],axis=0,ignore_index=True)

        else:
            print(f"Error in '{__name__}': Transcript dataframe invalid.")
            sys.exit(1)  # Exit the program with an error code (1)

    def parse(self,region=None,transcripts=None):

        if transcripts:
            self.transcripts = self.__process_transcript_file__(transcripts)

        self.vep_fields = self.vcf.header.info['CSQ'].description.split(":")[1].strip(" ").split("|")

        # get variants
        varDf = pd.DataFrame(columns=self.vep_fields + 'type filters chrom pos ref alt coverage altreads vaf'.split(' ')) 
        for variant in self.vcf.fetch(region=region):

            # get VEP annotation field
            vepCsq = variant.info.get('CSQ')
            if vepCsq is None:
                sys.exit("No VEP fields")

            df = pd.DataFrame(columns=self.vep_fields)
            for i in vepCsq:
                df = pd.concat([df,pd.DataFrame([dict(zip(df.columns,i.split("|")))])],axis=0,ignore_index=True)

            # if no symbol, use Gene ID
            df.loc[df['SYMBOL']=='','SYMBOL'] = df.loc[df['SYMBOL']=='','Gene']

            if len(variant.ref) == len(variant.alts[0]):
                df['type'] = 'SNV'
            else:
                df['type'] = 'INDEL'

            df['filters'] = ';'.join(variant.filter.keys())
            df['chrom'] = variant.chrom
            df['pos'] = variant.pos
            df['ref'] = variant.ref
            df['alt'] = variant.alts[0]
            df['coverage'] = variant.samples[0]['DP']
            df['altreads'] = variant.samples[0]['AD'][1]
            df['vaf'] = round(variant.samples[0]['AF'][0] * 100,4)

            varDf = pd.concat([varDf,df],axis=0,ignore_index=True)

        # merge in transcript info
        varDf = varDf.merge(self.transcripts,left_on='Feature',right_on='Transcript',how='left')

        # only variants w.r.t defined transcripts are reported.
        # if no variants are in these transcripts.
        varDf = varDf[~varDf['Transcript'].isna()].reset_index().drop(columns='index')
        if len(varDf) > 0:

            # get first consequence
            if 'Consequence' in varDf.columns:
                varDf['Consequence'] = varDf['Consequence'].apply(lambda x: x.split("&")[0])

            # intify distance
            if 'DISTANCE' in varDf.columns:
                varDf['DISTANCE'] = varDf['DISTANCE'].apply(lambda x: 0 if x=='' else int(x))

            if 'PICK' in varDf.columns:
                varDf['PICK'] = varDf['PICK'].apply(lambda x: 0 if x=='' else 1)

            if 'MAX_AF' in varDf.columns:
                varDf['MAX_AF'] = varDf['MAX_AF'].apply(lambda x: 0 if x=='' else round(float(x)*100,5))

            if 'HGVSc' in varDf.columns:
                varDf['csyntax'] = varDf['HGVSc'].apply(lambda x: x.split(':')[1] if ':' in x else 'noncoding')
                varDf.drop(columns='HGVSc')

                # *  ---->                        
                varDf.loc[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('upstream')) & (varDf['STRAND']==1),'csyntax'] = varDf[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('upstream')) & (varDf['STRAND']==1)].apply(lambda r: "c.-"+str(r['cdsStart'] - r['pos'])+r['ref']+'>'+r['alt'],axis=1)

                # <---- *
                varDf.loc[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('upstream')) & (varDf['STRAND']==-1),'csyntax'] = varDf[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('upstream')) & (varDf['STRAND']==1)].apply(lambda r: "c.-"+str(r['pos']-r['cdsEnd'])+revcomp(r['ref'])+'>'+revcomp(r['alt']),axis=1)

                # ---->  *                 
                varDf.loc[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('downstream')) & (varDf['STRAND']==-1),'csyntax'] = varDf[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('downstream')) & (varDf['STRAND']==1)].apply(lambda r: "c.+"+str(r['DISTANCE'])+r['ref']+'>'+r['alt'],axis=1)
        
                # *  <----
                varDf.loc[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('downstream')) & (varDf['STRAND']==-1),'csyntax'] = varDf[(varDf['csyntax']=='noncoding') & (varDf['Consequence'].str.contains('downstream')) & (df['STRAND']==1)].apply(lambda r: "c.+"+str(r['DISTANCE'])+r['ref']+'>'+r['alt'],axis=1)

            if 'HGVSp' in varDf.columns:
                varDf['psyntax'] = varDf['HGVSp'].apply(lambda x: re.sub("\%3D",'=',convert_aa(x.split(':')[1])) if ':' in x else None)
                varDf.loc[(varDf['psyntax'].isna()),'psyntax'] = varDf.loc[(varDf['psyntax'].isna()),'csyntax']
                varDf.drop(columns='HGVSp')

        varDf.rename(columns={'Gene':'GeneID','SYMBOL':'Gene'},inplace=True)
        self.df = varDf


class CnvVcf:
    def __init__(self,file=None):
        # initialize vcf object
        self.vcf = pysam.VariantFile(file,"r") if isinstance(file, str) else None
        self.__process_transcript_file__(transcripts)

        # initialize dataframe with small variants
        self.variants = pd.DataFrame(columns=['type','filters','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf'])

        self.df = pd.DataFrame()
        self.vep_fields = None

        if self.vcf != None:
            self.parse()


    def parse(self,file=None):

        cnvvcf = VCF(cnvvcffile)

        for variant in cnvvcf:
                
            vartype = variant.ALT
            if len(vartype) > 1:
                vartype = 'CNLOH'
            elif vartype[0] == '<DEL>':
                vartype = 'DEL'
            elif vartype[0] == '<DUP>':
                vartype = 'DUP'
            else:
                vartype = 'UNKNOWN'

            filter = 'PASS'
            if variant.FILTER is not None:
                filter = variant.FILTER

            chr1 = str(variant.CHROM)
            pos1 = variant.POS
            pos2 = variant.INFO.get('END')

            chr2 = chr1
            svlen = pos2 - pos1 + 1

            # get cytobands
            bands = 'None'
            bandstring = 'None'
            if variant.INFO.get('CYTOBANDS') is not None:
                bands = variant.INFO['CYTOBANDS'].split('&')
                bandstring = bands[0]
                if len(bands) > 1:
                    bandstring = bands[0] + bands[-1]

            # gene by overlap between variant and genes covqcdf. This is all we need for this resolution.
            genestring = 'None'
            genes = pr.PyRanges(covDf).intersect(pr.PyRanges(chromosomes = str(variant.CHROM),starts = [variant.POS],ends = [variant.INFO.get('END')])).df
            if genes.shape[0] > 0:
                genes = genes['Gene'].unique().tolist()
                if len(genes) == 0:
                    genestring = 'None'
                elif len(genes) > 10:
                    genestring = str(len(genes)) + " genes"
                else:
                    genestring = ",".join(genes)
                    
            # For example:  seq[GRCh38] del(X)(q21.31q22.1)
            #          chrX:g.89555676_100352080del

            csyntax = '.'
            psyntax = '.'
            if vartype == 'DEL':
                csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "del"
                if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                    psyntax = "seq[GRCh38] -" + chr1.replace('chr','')
                    
                elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
                    psyntax = "seq[GRCh38] -" + chr1.replace('chr','')

                elif bands[0].find('p') > -1:
                    bands.reverse()
                    psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                    
                else:
                    psyntax = "seq[GRCh38] del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            elif vartype == 'DUP':
                csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "dup"
                if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                    psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

                elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
                    psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

                elif bands[0].find('p') > -1:
                    bands.reverse()
                    psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                    
                else:
                    psyntax = "seq[GRCh38] dup(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            elif vartype == 'CNLOH':
                csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "cnLOH"
                if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                    psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

                elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
                    psyntax = "seq[GRCh38] +" + chr1.replace('chr','')

                elif bands[0].find('p') > -1:
                    bands.reverse()
                    psyntax = "seq[GRCh38] cnLOH(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                    
                else:
                    psyntax = "seq[GRCh38] cnLOH(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"

            # abundance
            abundance = variant.format('CF')[0][0]
            copynumber = variant.format('CN')[0][0]



class AssayQC:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame()
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []


class CovQC:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame(dataframe)
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []

class SV:
    def __init__(self, file=None, df=None):
        """
        Initialize the CovQC class with a pandas DataFrame, and two lists: genes and targets.

        :param dataframe: pandas DataFrame to be used in the class.
        :param genes: List of genes.
        :param targets: List of targets.
        """
        self.covdf = dataframe if isinstance(dataframe, pd.DataFrame) else pd.DataFrame(dataframe)
        self.genes = genes if isinstance(genes, list) else []
        self.targets = targets if isinstance(targets, list) else []

