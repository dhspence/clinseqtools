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

class ClinSeqTools:
    def __init__(self):

        # initialize dataframe with small variants
        self.variant_files = []
        self.variant_vep_fields = {}
        self.variant_df = pd.DataFrame()
        self.variants = pd.DataFrame(columns=['category','type','filters','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf'])

        self.sv_files = []
        self.sv_vep_fields = {}
        self.sv_df = pd.DataFrame()
        self.svs = pd.DataFrame(columns=['category','type','filters','chrom1','pos1','chrom2','pos2','genes','csyntax','psyntax','length','id','abundance','info'])

        self.qc_files = []
        self.qc_df = pd.DataFrame(columns=['category','metric','value','qcrange','flag','qcmetric'])

        # this is for exon/gene coverage metrics
        self.coverage_files = []
        self.coverage_df = pd.DataFrame(columns=['category','gene','type','region','mean','covered1','covered2'])

        # df for haplotect loci
        self.haplotect_files = []
        self.haplotect_df = pd.DataFrame(columns=['chr','pos1','pos2','pos1_allele1','pos1_allele2','pos2_allele1','pos2_allele2','pop_counts','distance','total_count','sample_counts','contam_fraction'])

        self.biomarker_files = []
        self.biomarkers = pd.DataFrame(columns=['metric','value'])

    def __process_transcript_file__(self,transcripts):
        if isinstance(transcripts, pd.DataFrame) and set(['Gene','Transcript','cdsStart','cdsEnd','Strand']).issubset(transcripts.columns):
            transcripts = transcripts.drop(columns='Gene')
            self.transcripts = pd.concat([transcripts[(transcripts['Strand']==1)].groupby(['Transcript','Strand']).agg({'cdsStart': 'min', 'cdsEnd': 'max'}).reset_index(),transcripts[(transcripts['Strand']==-1)].groupby(['Transcript','Strand']).agg({'cdsStart': 'max', 'cdsEnd': 'min'}).reset_index()],axis=0,ignore_index=True)

        else:
            print(f"Error in '{__name__}': Transcript dataframe invalid.")
            sys.exit(1)  # Exit the program with an error code (1)


    def parse_variant_file(self,vcf=None,region=None,transcripts=None):
        self.variant_files.append(vcf)
        vcf = pysam.VariantFile(vcf,"r") if isinstance(vcf, str) else None
        self.__process_transcript_file__(transcripts)

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

            """
            # do some error checking
            if len(df[(df['Gene']) & (df['SYMBOL']!=df['Gene'])])>0 or len(df[df['STRAND']!=df['Strand']])>0:
                print(f"Error in '{__name__}': Transcript information not congruent.")
                sys.exit(1)  # Exit the program with an error code (1)
            """

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

