import pysam
import pandas as pd
import sys, os, re

"""
Helper functions
"""

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

def sort_chrompos(row,chrom='chrom',pos='pos'):
    # Extract the numeric part from the 'Chromosome' value
    if row[chrom].startswith('chr'):
        num_part = row[chrom][3:]
        if num_part.isdigit():
            chromosome_value = int(num_part)
        else:
            chromosome_value = float('inf')
    else:
        chromosome_value = float('inf')
    
    return (chromosome_value, row[pos])

def vepToTable(csq,vep_fields):

    vep_fields = vep_fields.description.split(":")[1].strip(" ").split("|")
    if csq is None:
        sys.exit("No CSQ field found!")
        
    df = pd.DataFrame(columns=vep_fields)
    for i in csq:
        df = pd.concat([df,pd.DataFrame([dict(zip(df.columns,i.split("|")))])],axis=0,ignore_index=True)

    # if no symbol, use Gene ID
    df.loc[df['SYMBOL']=='','SYMBOL'] = df.loc[df['SYMBOL']=='','Gene']

    df['gene'] = df['SYMBOL']
    df['csyntax'] = ''
    df['psyntax'] = ''

    # get first consequence
    if 'Consequence' in df.columns:
        df['consequence'] = df['Consequence'].apply(lambda x: x.split("&")[0])

    # copy and rename selected columns.
    if 'Feature' in df.columns:
        df['transcript'] = df['Feature'].apply(lambda x: x.split("&")[0])

    if 'EXON' in df.columns:
        df['exon'] = df['EXON']

    # intify distance and pick
    if 'DISTANCE' in df.columns:
        df['DISTANCE'] = df['DISTANCE'].apply(lambda x: 0 if x=='' else int(x))

    if 'PICK' in df.columns:
        df['PICK'] = df['PICK'].apply(lambda x: 0 if x=='' else 1)

    if 'STRAND' in df.columns:
        df['STRAND'] = df['STRAND'].astype(int)

    if 'MAX_AF' in df.columns:
        df['pop_af'] = df['MAX_AF'].apply(lambda x: 0 if x=='' else round(float(x)*100,5))

    df = df.replace('', pd.NA)

    df['annotations'] = [[] for _ in range(len(df))]
    if 'Existing_variation' in df.columns:
        df['annotations'] = df['Existing_variation'].apply(lambda x: [] if pd.isna(x) else x.split('&'))

    # get assay-specific custom annotations
    dbcolumn = next((col for col in df.columns if '_DB' in col), None)
    if dbcolumn:
        df['annotations'] = df.apply(lambda r: r['annotations'] + ['DB=' + str(r[dbcolumn])] if pd.notna(r['annotations']) and pd.notna(r[dbcolumn]) else r['annotations'], axis=1)
        
    blacklist = next((col for col in df.columns if '_BLACKLIST' in col), None)
    if blacklist and (df[blacklist]==1).any():
        df['annotations'] = df.apply(lambda r: r['annotations'] + ['BLACKLIST'] if pd.notna(r['annotations']) and pd.notna(r[blacklist]) else r['annotations'], axis=1)

    whitelist = next((col for col in df.columns if '_WHITELIST' in col), None)
    if whitelist and (df[whitelist]==1).any():
        df['annotations'] = df.apply(lambda r: r['annotations'] + ['WHITELIST'] if pd.notna(r['annotations']) and pd.notna(r[whitelist]) else r['annotations'], axis=1)

    df['annotations'] = df['annotations'].apply(lambda x: ';'.join(x))

    return df

def vepSvGeneEffect(row):
    if row['EXON']:
        return('(exon' + row['EXON'].split('/')[0] + ')')
    elif row['INTRON']:
        return('(intron' + row['INTRON'].split('/')[0] + ')')
    elif row['DISTANCE']:
        if 'upstream' in row['Consequence']:
            return('(' + str(row['DISTANCE']) + 'bp upstream)')
        elif 'downstream' in row['Consequence']:
            return('(' + str(row['DISTANCE']) + 'bp downstream)')
        else:
            return('')
    else:
        return('')
        
class ClinSeqTools:
    def __init__(self):

        self.__version__ = '1.0.0'

        # set these static variables
        self.__acrocentrics__ = ['chr13','chr14','chr15','chr21','chr22']
        self.__consequences__ = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_ablation","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant"]

        # Small variant DF with all VEP fields.
        self.variant_fields = ['category','type','filters','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','pop_af','annotations','coverage','altreads','vaf']        
        # Standardized DF for small variants. Meant to be uniform across assays
        self.variants = pd.DataFrame(columns=self.variant_fields)
        self.variant_files = []

        # Standardized SV DF. Note: Gateway seq has a separate CNV table. We could  combine the CNV + SV or keep them separate.
        self.sv_fields = ['category','type','filters','chrom1','pos1','chrom2','pos2','genes','csyntax','psyntax','length','abundance','info']
        # Comprehensive Sv DF with all VEP fields.
        self.svs = pd.DataFrame(columns=self.sv_fields)
        self.sv_files = []

        # GatewaySeq CNV DF has this format: ['type','gene','chrom','start','end','copy_ratio','qual','filter','copynumber','ploidy','dblookup'])

        # Standardized QC DF. These come mainly from dragen. flag and qcmetric indicate whether the metric is out of range and whether the qc metric is tracked
        self.qc_df = pd.DataFrame(columns=['category','metric','value','qcrange','flag','qcmetric'])
        self.qc_files = []
        
        # Standardized DF for region coverage metrics. 'covered1' and 'covered2' indicate % of region at two coverage levels.
        self.coverage_df = pd.DataFrame(columns=['category','gene','type','region','mean','covered1','covered2'])
        self.coverage_files = []

        # Standardized DF for haplotect loci
        self.haplotect_df = pd.DataFrame(columns=['chr','pos1','pos2','pos1_allele1','pos1_allele2','pos2_allele1','pos2_allele2','pop_counts','distance','total_count','sample_counts','contam_fraction'])
        self.haplotect_files = []

        # Standardized DF for haplotect loci
        self.biomarkers = pd.DataFrame(columns=['metric','value'])
        self.biomarker_files = []

        self.transcripts = pd.DataFrame(columns=['gene','transcript','strand','cdsstart','cdsend'])


    """
    Internal functions
    """

    def _process_transcript_file(self,transcripts):
        if isinstance(transcripts, pd.DataFrame):
            if set(['gene','transcript','cdsstart','cdsend','strand']).issubset(transcripts.columns):
                transcripts = transcripts[transcripts['transcript']!='.'][['gene','transcript','cdsstart','cdsend','strand']]
                self.transcripts = pd.concat([transcripts[(transcripts['strand']=='1')].groupby(['gene','transcript','strand']).agg({'cdsstart': 'min', 'cdsend': 'max'}).reset_index(),transcripts[(transcripts['strand']=='-1')].groupby(['gene','transcript','strand']).agg({'cdsstart': 'max', 'cdsend': 'min'}).reset_index()],axis=0,ignore_index=True)

            else:
                print(f"Error in '{__name__}': Transcript dataframe is invalid.")
                sys.exit(1)  # Exit the program with an error code (1)

    def _vep_parse_sv_genes(self,variant,csq_header=None):
        # set up variant info DF
    
        vepDf = vepToTable(variant.info.get("CSQ"),csq_header)

        # annotate with known transcripts/genes, if given
        vepDf['KnownTrx'] = 0
        vepDf['KnownGene'] = 0
        if len(self.transcripts)>0:
            vepDf.loc[vepDf['Feature'].isin(self.transcripts['transcript'].tolist()),'KnownTrx'] = 1
            vepDf.loc[vepDf['SYMBOL'].isin(self.transcripts['gene'].tolist()),'KnownGene'] = 1

        # sort by known gene, distance, then transcript, then pick
        vepDf = vepDf.sort_values(by=['KnownGene','DISTANCE','KnownTrx','PICK'], ascending=[False,True,False,False]).drop_duplicates(subset='SYMBOL',keep='first')
        # annotates impact on gene (e.g,. full deletion or partial)
        vepDf['GeneEffect'] = vepDf.apply(lambda r: vepSvGeneEffect(r),axis=1)
        # determines whether consequence is significant (e.g., deletes a gene vs. is upstream)
        vepDf['GeneImpact'] = vepDf['Consequence'].apply(lambda r: int(len(set(r.split('&')) & set(self.__consequences__))>0))

        # After all this annotation, we only get the events that affects a known gene
        vepDf = vepDf[vepDf['KnownGene']==1] #[['SYMBOL','GeneImpact','GeneEffect']]
        if 'KnownSvGenes' in variant.info.keys():
            if variant.info.get('KnownSvGenes'):
                knownGeneDf = pd.concat([knownGeneDf,pd.DataFrame({'SYMBOL':variant.info.get('KnownSvGenes').split(','),'GeneEffect':'','GeneImpact':1})])

        vepDf.sort_values(by=['SYMBOL','GeneImpact','GeneEffect'],ascending=[True,False,True],na_position='last').drop_duplicates(subset='SYMBOL',keep='first')

        vepDf['category'] = ''
        vepDf['type'] = variant.info.get('SVTYPE')
        vepDf['filters'] = ''
        vepDf['chrom1'] = variant.chrom
        vepDf['pos1'] = variant.pos
        vepDf['chrom2'] = variant.chrom
        vepDf['pos2'] = variant.info.get('END')

        if len(variant.filter.keys())==0:
            vepDf['filters'] = 'PASS'
        else:
            vepDf['filters'] = ';'.join(variant.filter.keys())

        vepDf['length'] = None
        if 'SVLEN' in variant.info.keys():
            vepDf['length'] = variant.info.get('SVLEN')

        return vepDf
    
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
        self.coverage_df = pd.concat([covDf.iloc[:,:4],ids,covDf.iloc[:,5:]],axis=1)
        self.coverage_df['cdsstart'] = self.coverage_df['cdsstart'].astype(int)
        self.coverage_df['cdsend'] = self.coverage_df['cdsend'].astype(int)
        self._process_transcript_file(self.coverage_df)
        self.coverage_files.append(file)

    def parse_variant_file(self,file=None,transcripts=None,extra_fields=None):       

        # add these vep/info fields to the standard list 
        if extra_fields:
            self.variant_fields.append(extra_fields)
            
        self.variant_files.append(checkfile(file))
        vcf = pysam.VariantFile(file,"r")

        if isinstance(transcripts,pd.DataFrame):
            self._process_transcript_file(transcripts)

        # get variants
        varDf = pd.DataFrame(columns=self.variant_fields)
        for variant in vcf.fetch():

            # get vep annotations
            df = vepToTable(variant.info.get("CSQ"),vcf.header.info['CSQ'])
            # add info annotations
            for f in variant.header.info.keys():
                # note, this will replace annotations in the VEP header
                if isinstance(variant.info.get(f), (list, tuple)):
                    df[f] = ','.join(map(str, variant.info.get(f)))
                else:
                    df[f] = variant.info.get(f)
                
            # add variant info
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

            df['cds_position'] = df['DISTANCE']

            # merge annotation with assay transcripts, if given.
            # Then find relative cds position (for upstream/downstream variants)
            if len(self.transcripts)>0:
                df = df.merge(self.transcripts.drop(columns='gene'),on='transcript',how='inner')
                
                df.loc[(df['consequence'].str.contains('upstream')) & (df['STRAND']==1),'cds_position'] = df[((df['consequence'].str.contains('upstream')) & (df['STRAND']==1))].apply(lambda r: r['cdsstart']-r['pos'],axis=1)
                df.loc[(df['consequence'].str.contains('upstream')) & (df['STRAND']==-1),'cds_position'] = df[((df['consequence'].str.contains('upstream')) & (df['STRAND']==-1))].apply(lambda r: r['pos']-r['cdsend'],axis=1)

                df.loc[(df['consequence'].str.contains('downstream')) & (df['STRAND']==1),'cds_position'] = df[((df['consequence'].str.contains('downstream')) & (df['STRAND']==1))].apply(lambda r: r['pos']-r['cdsend'],axis=1)
                df.loc[(df['consequence'].str.contains('downstream')) & (df['STRAND']==-1),'cds_position'] = df[((df['consequence'].str.contains('downstream')) & (df['STRAND']==-1))].apply(lambda r: r['cdsend']-r['pos'],axis=1)

            # create c and p syntaxes
            if 'HGVSc' in df.columns:
                df['csyntax'] = df['HGVSc'].apply(lambda x: x.split(':')[1] if pd.notna(x) and ':' in x else 'noncoding')
                df.drop(columns='HGVSc',inplace=True)

                # create valid upstream/downstream csyntax
                # *  ---->                        
                df.loc[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('upstream')) & (df['STRAND']==1),'csyntax'] = df[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('upstream')) & (df['STRAND']==1)].apply(lambda r: "c.-"+str(r['cds_position'])+r['ref']+'>'+r['alt'],axis=1)

                # <---- *
                df.loc[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('upstream')) & (df['STRAND']==-1),'csyntax'] = df[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('upstream')) & (df['STRAND']==-1)].apply(lambda r: "c.-"+str(r['cds_position'])+revcomp(r['ref'])+'>'+revcomp(r['alt']),axis=1)

                # ---->  *                 
                df.loc[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('downstream')) & (df['STRAND']==1),'csyntax'] = df[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('downstream')) & (df['STRAND']==1)].apply(lambda r: "c.+"+str(r['cds_position'])+r['ref']+'>'+r['alt'],axis=1)
        
                # *  <----
                df.loc[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('downstream')) & (df['STRAND']==-1),'csyntax'] = df[(df['csyntax']=='noncoding') & (df['consequence'].str.contains('downstream')) & (df['STRAND']==-1)].apply(lambda r: "c.+"+str(r['cds_position'])+revcomp(r['ref'])+'>'+revcomp(r['alt']),axis=1)

            if 'HGVSp' in df.columns:
                df['psyntax'] = df['HGVSp'].apply(lambda x: re.sub("\%3D",'=',convert_aa(x.split(':')[1])) if pd.notna(x) and ':' in x else None)
                df.loc[(df['psyntax'].isna()),'psyntax'] = df.loc[(df['psyntax'].isna()),'csyntax']
                df.drop(columns='HGVSp',inplace=True)

            df = df.sort_values(by=['chrom','pos','ref','alt','transcript','PICK'],ascending=[True,True,True,True,True,False],na_position='last').drop_duplicates(subset=['chrom','pos','ref','alt'],keep='first')

            # sort variants and keep those associated with a known trx or have the pick flag
            varDf = pd.concat([varDf,df[self.variant_fields]],axis=0,ignore_index=True)

        self.variants = varDf
    
    def get_filtered_variants(self,consequences=None,pop_af=0.01,filter_condition=None):
        if isinstance(consequences,list):
            self.__consequences__ = consequences

        if isinstance(pop_af,float):
            self.__pop_af__ = pop_af
        
        condition = (((self.variants['pop_af']<self.__pop_af__) & 
                     (self.variants['consequence'].isin(self.__consequences__)) &
                     (~self.variants['annotations'].str.contains('BLACKLIST',na=False))) |
                     (self.variants['annotations'].str.contains('WHITELIST',na=False)))
        
        return self.variants[condition]
    """
    def parse_sv_file(self,file=None,recurrentsvs=None,transcripts=None,extra_fields=None):

        # add extra columns 
        if isinstance(extra_fields,str):
            extra_fields = [extra_fields]

        # get known fusion events
        if recurrentsvs and checkfile(recurrentsvs):
            recurrentsvs = pd.read_csv(recurrentsvs,sep=',', header=None)            
            recurrentsvs.columns = ['gene1','strand1','gene2','strand2']

        self.sv_files.append(checkfile(file))
        svvcf = pysam.VariantFile(file,"r")

        if isinstance(transcripts,pd.DataFrame):
            self._process_transcript_file(transcripts)
            
        # go through all variants once and store them before parsing.
        # this is because we need both ends of the BNDs

        vardict = {}
        alreadydone = set() # set holding completed variants (for BNDs that have 2 records)

        for variant in svvcf.fetch():
            
            # save all BNDs because we dont know if one end
            # hits a recurrent gene
            if variant.info.get('SVTYPE')=='BND':
                vardict[variant.id] = variant

            else:
                # if not a BND, store variants that PASS or invovle a known gene
                if len(variant.filter.keys())==0 and 'KnownSvGenes' in variant.info.keys() and variant.info.get('KnownSvGene') is not None:
                    vardict[variant.id]

        # iterate over all variants to parse/format them
        varDf = pd.DataFrame(columns=self.sv_fields)

        for v in vardict.keys():

            variant = vardict[v]
            mate = None

            # get mate for BNDs
            if variant.info['SVTYPE'] == 'BND':
                if variant.info.get('MATEID') in alreadydone or not variant.info.get('MATEID') in vardict:
                    continue
                else:
                    mate = vardict[variant.info.get('MATEID')]
                    # if KnownSvGene is not present and the filter is not PASS for both variant and mate, then skip
                    if (len(variant.filter.keys())>0 and
                        'KnownSvGenes' in variant.info.keys() and
                            not variant.info.get('KnownSvGene') and
                                len(mate.filter.keys())>0 and
                        'KnownSvGenes' in mate.info.keys() and
                            not mate.info.get('KnownSvGene')):
                        continue
                
                # this creates a comprehensive DF with all variant info, including vep fields.
                # it only retains one entry per gene, selected based on assay info.
                varDf = self._vep_parse_sv_genes(variant,csq_header=svvcf.header.info['CSQ'])
                
                # Format DEL/DUP/INS
                if variant.info['SVTYPE'] in ['DEL','DUP','INS']:
                    cytobands = varDf['cytobands'][0].split("&")
                
                    # Determine if this is a DEL/DUP that is recurrent
                    # -Criteria: 1 or 2 known genes affected by this event and they match events in our list
                    # -Since these are del/dup/ins events, the strands must be the same
                    geneHits = varDf[varDf['GeneImpact']==1]['SYMBOL'].tolist()
                    genePairHit = recurrentsvs[recurrentsvs.gene1.isin(geneHits) & recurrentsvs.gene2.isin(geneHits)]
                    
                    if genePairHit.shape[0] > 0 and genePairHit[genePairHit['strand1']==genePairHit['strand2']].shape[0] > 0:
                        varDf['category'] = 'RECURRENTSV'
                    elif len(variant.filters.keys())==0:
                        varDf['category'] = 'OTHERSV'
                
                    var = varDf.iloc[0].to_dict()
                    # format c and p syntax
                    if var['type'] == 'DEL':
                        var['csyntax'] = f"{var['chrom1']}:g.{var['pos1']}_{var['pos2']}del"

                        if cytobands[0].find('p') > -1 and cytobands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                            var['psyntax'] = "seq[GRCh38] -" + var['chrom1'].replace('chr','')
                            
                        elif 'q11' in cytobands and 'qter' in cytobands and var['chrom1'] in self.__acrocentrics__:
                            var['psyntax'] = "seq[GRCh38] -" + var['chrom1'].replace('chr','')

                        elif cytobands[0].find('p') > -1:
                            cytobands.reverse()
                            var['psyntax'] = "seq[GRCh38] del(" + var['chrom1'].replace('chr','') + ")(" + cytobands[0] + cytobands[-1] + ")"
                            
                        else:
                            var['psyntax'] = "seq[GRCh38] del(" + var['chrom1'].replace('chr','') + ")(" + cytobands[0] + cytobands[-1] + ")"
                    
                    elif var['type'] == 'DUP' or var['type']== 'INS':
                        var['csyntax'] = var['chrom1'] + ":g." + str(var['pos1']) + "_" + str(var['pos2']) + var['type'].lower()

                        if cytobands[0].find('p') > -1 and cytobands[-1].find('q') > -1:
                            var['psyntax'] = "seq[GRCh38] +" + var['chrom1'].replace('chr','')

                        elif 'q11' in cytobands and 'qter' in cytobands and varDf['chrom1'] in self.__acrocentrics__:
                            var['psyntax'] = "seq[GRCh38] +" + var['chrom1'].replace('chr','')

                        elif cytobands[0].find('p') > -1:
                            cytobands.reverse()
                            var['psyntax'] = "seq[GRCh38] " + var['type'].lower() + "(" + var['chrom1'].replace('chr','') + ")(" + cytobands[0] + cytobands[-1] + ")"
                            
                        else:
                            var['psyntax'] = "seq[GRCh38] " + var['type'].lower() + "(" + var['chrom1'].replace('chr','') + ")(" + cytobands[0] + cytobands[-1] + ")"

                SR = (0,0)
                PR = (0,0)
                if 'PR' in variant.format.keys():
                    PR = variant.samples[0]['PR']

                if 'SR' in variant.format.keys():
                    SR = variant.samples[0]['SR']

                varDf['abundance'] = (SR[1] + PR[1]) / (PR[0] + PR[1] + SR[0] + SR[1])

                varDf['info'] = ['ID=' + variant.id, 'PR=' + str(PR[1]) + '/' + str(PR[0]+PR[1]), 'SR=' + str(SR[1]) + '/' + str(SR[0]+SR[1]),'CONTIG=' + str(variant.info.get('CONTIG'))]

                for f in extra_fields:
                    if f in knownGeneDf.columns:
                        varDf['info'].append(f"{f}={','.join(knownGeneDf[f].tolist())}")
                    elif f in variant.info.keys() and variant.info.get(f):
                        if isinstance(variant.info.get(f),(list,tuple)):
                            varDf['info'].append(f"{f}={','.join(variant.info.get(f))}")
                        else:
                            varDf['info'].append(f"{f}={variant.info.get(f)}")

                varDf['info'] = varDf['info'].apply(lambda x: ';'.join(varDf['info']))

                self.svs = pd.concat([self.svs,varDf[[self.svs.columns]]],axis=0,ignore_index=True)

        
        # now handle BNDs, which each have 2 entries.
        # this includes translocations and inversions
        alreadydone = set()
        i  = 0
        for v in passedvars.items():

            # get the first record
            variant = v[1]

            # skip if this is the mate and we already processed the record pair
            if variant.info.get('MATEID') in alreadydone or variant.info.get('MATEID') not in passedvars:
                continue

            # get the mate
            mate = passedvars[variant.info.get('MATEID')]

            # initial filtering, to speed things up:
            # if KnownSvGene is not present and the filter is not PASS, then skip
            if len(variant.filter.keys())>0 and len(mate.filter.keys())>0:
                if 'KnownSvGenes' in variant.info.keys() and not variant.info.get('KnownSvGene') and not mate.info.get('KnownSvGene'):
                    continue

            varDf = pd.DataFrame(columns=self.svs.columns)

            varDf['category'] = ''
            varDf['type'] = variant.info.get('SVTYPE')

            # check to see if this is an inversion
            if variant.info.get('INV3') is not None or variant.info.get('INV5') is not None:
                varDf['type'] = 'INV'
            
            # get filter info
            if len(variant.filter.keys())>0:
                varDf['filters'] = ';'.join(variant.filter.keys())
            else:
                varDf['filters'] = 'PASS'
            
            # get gene 1 info
            varDf['chrom1'] = str(variant.chrom)
            varDf['pos1'] = variant.pos
            varDf['chrom2'] = mate.chrom
            varDf['pos2'] = mate.pos

            genes1='NA'
            transcript1='NA'
            region1 = 'NA'
            strand1 = '+'
            bands1 = 'NA'

            genes2='NA'
            transcript2='NA'
            region2 = 'NA'
            strand2 = '+'
            bands2 = 'NA'

            # get gene info for VARIANT
            vepDf = vepToTable(variant.info.get("CSQ"),svvcf.header.info['CSQ'])

            vepDf['KnownTrx'] = 0
            vepDf['KnownGene'] = 0
            if len(self.transcripts) > 0:
                vepDf.loc[vepDf['Feature'].isin(self.transcripts['transcript'].tolist()),'KnownTrx'] = 1
                vepDf.loc[vepDf['SYMBOL'].isin(self.transcripts['gene'].tolist()),'KnownGene'] = 1

            bands1 = vepDf['cytobands'][0].split("&")[0]

            vepCsq = vepCsq.sort_values(by=['KnownGene','DISTANCE','KnownTrx','PICK'], ascending=[False,True,False,False]).drop_duplicates(subset='SYMBOL',keep='first')
            vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepSvGeneEffect(r),axis=1)
            gene1Df= vepCsq[['SYMBOL','KnownGene','KnownTrx','DISTANCE','STRAND','GeneEffect']]

            if variant.INFO.get('KnownSvGenes') is not None:
                gene1Df = pd.concat([gene1Df,pd.merge(pd.DataFrame({'SYMBOL':variant.INFO['KnownSvGenes'].split(','),'KnownGene':1}),
                                pd.concat([recurrentsvs[['gene1','strand1']].rename(columns={'gene1':'SYMBOL','strand1':'STRAND'}),recurrentsvs[['gene2','strand2']].rename(columns={'gene2':'SYMBOL','strand2':'STRAND'})],axis=0).drop_duplicates(),on='SYMBOL')],axis=0)

            # sort by known transcript and then distance and get only the first transcript.
            gene1Df = gene1Df.sort_values(by=['KnownGene','DISTANCE','KnownTrx','GeneEffect'], ascending=[False,True,False,False],na_position='last').head(1).fillna('')
            
            if gene1Df.shape[0] > 0:
                genes1 = gene1Df['SYMBOL'].tolist()[0]
                region1 = gene1Df['GeneEffect'].tolist()[0]
                strand1 = gene1Df['STRAND'].tolist()[0]

            # get gene info for MATE
            vepCsq = vepToTable(mate.INFO['CSQ'],svvcf.get_header_type('CSQ'))    
            vepCsq['KnownTrx'] = 0
            vepCsq.loc[vepCsq['Feature'].isin(knownTrx['Transcript'].tolist()),'KnownTrx'] = 1
            vepCsq['KnownGene'] = 0
            vepCsq.loc[vepCsq['SYMBOL'].isin(knownTrx['Gene'].tolist()),'KnownGene'] = 1

            bands2 = vepCsq['cytobands'][0].split("&")[0]

            vepCsq = vepCsq.sort_values(by=['KnownGene','DISTANCE','KnownTrx','PICK'], ascending=[False,True,False,False]).drop_duplicates(subset='SYMBOL',keep='first')
            vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
            gene2Df= vepCsq[['SYMBOL','KnownGene','KnownTrx','DISTANCE','STRAND','GeneEffect']]

            if mate.INFO.get('KnownSvGenes') is not None:
                gene2Df = pd.concat([gene2Df,pd.merge(pd.DataFrame({'SYMBOL':mate.INFO['KnownSvGenes'].split(','),'KnownGene':1}),
                                pd.concat([recurrentSvs[['gene1','strand1']].rename(columns={'gene1':'SYMBOL','strand1':'STRAND'}),recurrentSvs[['gene2','strand2']].rename(columns={'gene2':'SYMBOL','strand2':'STRAND'})],axis=0).drop_duplicates(),on='SYMBOL')],axis=0)

            # sort by known transcript and then distance and get only the first transcript.
            gene2Df = gene2Df.sort_values(by=['KnownGene','DISTANCE','KnownTrx','GeneEffect'], ascending=[False,True,False,False],na_position='last').head(1).fillna('')

            if gene2Df.shape[0] > 0:
                genes2 = gene2Df['SYMBOL'].tolist()[0]
                region2 = gene2Df['GeneEffect'].tolist()[0]
                strand2 = gene2Df['STRAND'].tolist()[0]

            if genes1=='NA' or genes1=='':
                genes1 = 'INTERGENIC'

            if genes2=='NA' or genes2=='':
                genes2 = 'INTERGENIC'

            orientation = '+'
            if variant.ALT[0].find("[") == 0 or variant.ALT[0].find("]") > 0:
                orientation = '-'    

            # abundance
            abundance = 0.0
            PR = (0,0)
            SR = (0,0)            
            if variant.format("SR") is not None:
                SR = variant.format("SR")[0]
                    
            if variant.format("PR")[0] is not None:                
                PR =  variant.format("PR")[0]

            abundance = round((SR[1] + PR[1]) / (PR[0] + PR[1] + SR[0] + SR[1]) * 100,1)

            alt = variant.ALT[0]
            orientation = '+'
            if alt.find("[") == 0 or alt.find("]") > 0:
                orientation = '-'

            csyntax = 'NA'
            psyntax = 'NA'
            genestring = 'NA'

            if (chr1.find('M') == -1 and chr2.find('M') == -1 and chr1.find('X') == -1 and chr2.find('X') == -1 and chr1.find('Y') == -1 and chr2.find('Y') == -1 and int(chr1.replace('chr','')) < int(chr2.replace('chr',''))) or chr1.find('X') > -1 or chr1.find('Y') > -1:
                csyntax = chr1 + ":g." + str(pos1) + "(+)::" + chr2 + ":g." + str(pos2) + "(" + orientation + ")"
                psyntax = 'seq[GRCh38] t(' + chr1.replace('chr','') + ';' + chr2.replace('chr','') + ')(' + bands1 + ';' + bands2 + ')'
                genestring = genes1+'('+strand1+')'+region1+'::'+genes2+'('+strand2+')'+region2
            else:
                csyntax = chr2 + ":g." + str(pos2) + "(+)::" + chr1 + ":g." + str(pos1) + "(" + orientation + ")"
                psyntax = 'seq[GRCh38] t(' + chr2.replace('chr','') + ';' + chr1.replace('chr','') + ')(' + bands2 + ';' + bands1 + ')'
                genestring = genes2+'('+strand2+')'+region2+'::'+genes1+'('+strand1+')'+region1       

            infostring = 'PR_READS=' + str(PR[1]) + '/' + str(PR[0]+PR[1]) + ';SR_READS=' + str(SR[1]) + '/' + str(SR[0]+SR[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))

            if len(filter) == 0:
                filter = 'PASS'
            else:
                filter = ';'.join(filter)

            isRecurrentSv = False

            # if one side of the BND involves a gene that can partner with any other gene and its a PASS variant then call it recurrent
            if recurrentSvs[(recurrentSvs['gene1']==genes1) & (recurrentSvs['gene2']=='*')].shape[0]>0 and genes1!=genes2 and genes2!='INTERGENIC' and filter=='PASS':
                isRecurrentSv = True

            if recurrentSvs[(recurrentSvs['gene1']==genes2) & (recurrentSvs['gene2']=='*')].shape[0]>0 and genes1!=genes2 and genes1!='INTERGENIC' and filter=='PASS':
                isRecurrentSv = True

            # if the genes and orientation match and the orientation of at least one gene doesnt matter (like in IgH rearrangements)
            if recurrentSvs[(recurrentSvs['gene1']==genes1) & (recurrentSvs['gene2']==genes2) & ((recurrentSvs['strand1']=='*') | (recurrentSvs['strand2']=='*'))].shape[0]>0:
                isRecurrentSv = True

            # else if the genes and orientation match a recurrent event
            if orientation == '+' and recurrentSvs[(recurrentSvs['gene1']==genes1) & (recurrentSvs['gene2']==genes2) & (recurrentSvs['strand1']==recurrentSvs['strand2'])].shape[0]>0:
                isRecurrentSv = True

            if orientation == '+' and recurrentSvs[(recurrentSvs['gene1']==genes2) & (recurrentSvs['gene2']==genes1) & (recurrentSvs['strand1']==recurrentSvs['strand2'])].shape[0]>0:
                isRecurrentSv = True

            if orientation == '-' and recurrentSvs[(recurrentSvs['gene1']==genes1) & (recurrentSvs['gene2']==genes2) & (recurrentSvs['strand1']!=recurrentSvs['strand2'])].shape[0]>0:
                isRecurrentSv = True

            if orientation == '-' and recurrentSvs[(recurrentSvs['gene1']==genes2) & (recurrentSvs['gene2']==genes1) & (recurrentSvs['strand1']!=recurrentSvs['strand2'])].shape[0]>0:
                isRecurrentSv = True

            # skip non-PASS variants unless they are recurrent
            if filter != 'PASS' and isRecurrentSv is False:
                continue

            # categorize and filter variants
            category = ''
            
            # report all recurrent events
            if isRecurrentSv is True: 
                category = 'Recurrent'

            # PASS variants involving known genes
            elif filter=='PASS' and knownTrx['Gene'].isin([genes1,genes2]).any():
                category = 'OtherSv'

            # PASS variants where both ends involve genes
            elif filter=='PASS' and genes1!='INTERGENIC' and genes2!='INTERGENIC':
                category = 'OtherSv'

            else:
                continue

            svs = pd.concat([svs,pd.DataFrame([dict(zip(svs.columns,[category,vartype,chr1,str(pos1),chr2,str(pos2),str(svlen),csyntax,psyntax,genestring,filter,str(variant.ID) + ";" + str(mate.ID),str(abundance)+"%",infostring,'NA']))])])

            alreadydone.add(variant.ID)
            """