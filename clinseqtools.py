import pysam
import pandas as pd
import sys, os, re

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    return file_path

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

        self.__version__ = '1.0.0'

        self.__acrocentrics__ = ['chr13','chr14','chr15','chr21','chr22']

        # initialize dataframe with small variants
        self.variant_files = []
        self.variant_vep_fields = {}

        # this is a comprehensive small variant DF with all VEP fields.
        self.variants_df = pd.DataFrame()
        
        # this is a standardized DF for small variants. Meant to be uniform across assays
        self.variants = pd.DataFrame(columns=['category','type','filters','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','max_af','annotations','coverage','altreads','vaf'])

        self.sv_files = []
        self.sv_vep_fields = {}

        # Comprehensive Sv DF with all VEP fields.
        self.svs_df = pd.DataFrame()

        # Standardized SV DF. Note: Gateway seq has a separate CNV table. We could  combine the CNV + SV or keep them separate.
        self.svs = pd.DataFrame(columns=['category','type','filters','chrom1','pos1','chrom2','pos2','genes','csyntax','psyntax','length','id','abundance','info'])
                                        #['type','gene','chrom','start','end','copy_ratio','qual','filter','copynumber','ploidy','dblookup'])
        self.transcript_df = pd.DataFrame(columns=['gene','transcript','strand','cdsstart','cdsend'])

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

    """
    Internal functions
    """

    def process_transcript_file(self,transcripts):
        if isinstance(transcripts, pd.DataFrame):
            if set(['transcript','cdsstart','cdsend','strand']).issubset(transcripts.columns):
                if 'gene' in transcripts.columns:
                    transcripts = transcripts.drop(columns='gene')
                     
                self.transcripts = pd.concat([transcripts[(transcripts['strand']==1)].groupby(['transcript','strand']).agg({'cdsstart': 'min', 'cdsend': 'max'}).reset_index(),transcripts[(transcripts['strand']==-1)].groupby(['transcript','strand']).agg({'cdsstart': 'max', 'cdsend': 'min'}).reset_index()],axis=0,ignore_index=True)

            else:
                print(f"Error in '{__name__}': Transcript dataframe invalid.")
                sys.exit(1)  # Exit the program with an error code (1)


    """
    External functions
    """

    def parse_coverage_report(self,file=None):
        with open(checkfile(file), 'r') as f:
            lines = f.readlines()

        # Remove '#' from the first line and use it as the header
        header_line = lines[0].replace('#', '').strip()
        covDf = pd.read_csv(file, header=None, skiprows=1, names=['chromosome','start','end','gene','info'] + header_line.split('\t')[5:],sep="\t")

        ids = covDf['info'].str.split('|',expand=True).replace('\s\+\s\S+','',regex=True)
        ids.columns = ['group','geneid','transcript','region','cdsstart','cdsend','strand']
        self.coverage_files.append(file)
        self.coverage_df = pd.concat([covDf.loc[:,:5],ids,covDf.loc[:,5:]],axis=1)

    def parse_variant_file(self,file=None,region=None,transcripts=None):       
            
        self.variant_files.append(checkfile(file))
        vcf = pysam.VariantFile(file,"r")

        if isinstance(transcripts,pd.DataFrame):
            self.process_transcript_file(transcripts)

        self.variant_vep_fields = vcf.header.info['CSQ'].description.split(":")[1].strip(" ").split("|")

        # get variants
        varDf = pd.DataFrame(columns=self.variant_vep_fields + 'type filters chrom pos ref alt coverage altreads vaf'.split(' ')) 
        for variant in vcf.fetch(region=region):

            # get VEP annotation field
            vepCsq = variant.info.get('CSQ')
            if vepCsq is None:
                sys.exit("No VEP fields")

            df = pd.DataFrame(columns=self.variant_vep_fields)
            for i in vepCsq:
                df = pd.concat([df,pd.DataFrame([dict(zip(df.columns,i.split("|")))])],axis=0,ignore_index=True)

            # if no symbol, use Gene ID
            df.loc[df['SYMBOL']=='','SYMBOL'] = df.loc[df['SYMBOL']=='','Gene']

            df['category'] = ''

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
        varDf = varDf.merge(self.transcripts,left_on='Feature',right_on='transcript',how='left')

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
            varDf.drop(columns='HGVSc',inplace=True)

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
            varDf.drop(columns='HGVSp',inplace=True)

        # reformat df and return
        varDf.rename(columns={'Gene':'GeneID','SYMBOL':'gene'},inplace=True)
        varDf.drop(columns=['Allele','strand'],inplace=True)
        varDf.columns = varDf.columns.str.lower()

        # sort variants and keep those associated with a known trx or have the pick flag (then its arbitrary)
        varDf = varDf.sort_values(by=['chrom','pos','ref','alt','transcript','pick'],ascending=[True,True,True,True,True,False],na_position='last').drop_duplicates(subset=['chrom','pos','ref','alt'],keep='first')
        self.variants_df = varDf

