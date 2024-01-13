#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse, subprocess, gzip
import sqlite3
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path

__version__ = '2.0.0'

ACROCENTRICS = ['13','14','15','21','22']

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def fileexists(file_path):
    """Check if a file exists at the given path."""
    if os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The outfile {file_path} exists!")
    return file_path

def revcomp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}  # DNA complement pairs
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna))
    return reverse_complement

def pos2codon(exonstart,exonend,cpos,pos,strand):
    if pos<=exonend and pos>=exonstart:
        if strand == '+':
            return(int((pos-exonstart + cpos) / 3 + .99))

        elif strand=='-':
             return(int((exonend-pos + cpos) / 3 + .99))

def decode_hex(string):
    hex_string = string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

def convert_aa(codon):
    three = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"]
    one  = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    
    for i in range(0,len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i],codon)

    return codon

def sort_chrompos(row,chrom='Chromosome',pos='Start'):
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

def checkQc(value, minmax):
    if len(minmax)==1:
        if value < float(minmax[0]):
            return ['>' + str(minmax[0]),'(!)']
        else: 
            return ['>' + str(minmax[0]),'']
    elif len(minmax)==2:
        if value > float(minmax[0]) and value < float(minmax[1]):
            return [str(minmax[0])+'-'+str(minmax[1]),'']
        else:
            return [str(minmax[0])+'-'+str(minmax[1]), '(!)']
        
    else:
        sys.exit("Reference range must have minimum and or maximum")

def vepToTable(csq,header):
    fields = header['Description'].strip('"').split("|")
    df = pd.DataFrame(columns=fields)
    for i in csq.split(','):
        df = pd.concat([df,pd.DataFrame([dict(zip(df.columns,i.split("|")))])],axis=0,ignore_index=True)

    # if no symbol, use Gene ID
    df.loc[df['SYMBOL']=='','SYMBOL'] = df.loc[df['SYMBOL']=='','Gene']

    df.loc[df['STRAND']=='1','STRAND'] = '+'
    df.loc[df['STRAND']=='-1','STRAND'] = '-'

    if 'DISTANCE' in df.columns:
        df['DISTANCE'] = df['DISTANCE'].apply(lambda x: 0 if x=='' else int(x))

    if 'PICK' in df.columns:
        df['PICK'] = df['PICK'].apply(lambda x: 0 if x=='' else 1)

    return df

def getVepFields(vcf):
    vep = {}
    i = 0
    for j in vcf.get_header_type('CSQ')['Description'].strip('"').split("|"):
        vep[j] = i
        i+=1
    return(vep)

def vepGeneEffect(row):
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
    
def dbAnnotateVariants(db,queryPass,queryFiltered,params):                
    cursor = db.execute(queryPass, params)
    countPass = cursor.fetchone()[0]

    cursor = db.execute(queryFiltered, params)
    countFiltered = cursor.fetchone()[0]

    return 'PASS='+str(countPass)+'/'+'Filtered='+str(countFiltered)

def dbUpdateVariants(db,table,info,vardf):
     # first remove existing records for this case
    query = "DELETE FROM " + table + " WHERE name = ? AND mrn = ? AND accession = ?"
    cursor = db.cursor()
    cursor.execute(query, (info['name'],info['mrn'],info['accession']))
    db.commit()
    # now add records
    df = vardf.copy()
    for k in info.keys():
        df[k] = info[k]

    df.to_sql(table, db, if_exists='append', index=False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Script
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description='Make ChromoSeq report')
parser.add_argument('-n','--name',required=True,type=str,help='Sample name')
parser.add_argument('-d','--dir',required=True,type=str,help='Output directory')
parser.add_argument('--qc',required=True,help='Assay spec json file')
parser.add_argument('--genebed',type=checkfile,help='Gene bed file')
parser.add_argument('--svbed',type=checkfile,help='SV bed file')
parser.add_argument('--svtargets',type=checkfile,help='Recurrent SV target list')
parser.add_argument('-p','--popaf',default=0.1,type=float,help='Max population AF to report small variants (%)')
parser.add_argument('-m','--mrn',default='NONE',type=str,help='Sample MRN number')
parser.add_argument('-a','--accession',default='NONE',type=str,help='Sample accession number')
parser.add_argument('-s','--specimen',default='NONE',type=str,help='Sample specimen type')
parser.add_argument('-b','--DOB',default='NONE',help='Date of birth')
parser.add_argument('--sex',default='NONE',help='Case sex for chromosome analysis')
parser.add_argument('-e','--exception',default='NONE',help='Exception')
parser.add_argument('-r','--runinfostr',default='NONE',help='Illumina Run Information String')
parser.add_argument('-V','--assayversion',default='NONE',help='Assay version')
parser.add_argument('--variantdb',required=False,type=checkfile,default=None,help='Sqlite database file for variant lookup and recording')
parser.add_argument('--noupdate',required=False,default=False,action='store_true',help='Dont update variant database')
parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

args = parser.parse_args()

# get accessory files (as input or from repo location)
geneBed = args.genebed
if not Path(geneBed).is_file():
    sys.exit("Gene bed file " + str(geneBed) + " not valid.")

svBed = args.svbed
if not Path(svBed).is_file():
    sys.exit("SV bed file " + str(svBed) + " not valid.")

svGenes = args.svtargets
if not Path(svGenes).is_file():
    sys.exit("SV target gene list " + str(svGenes) + " not valid.")

qcFile = args.qc
if not Path(qcFile).is_file():
    sys.exit("QC Specification file " + str(qcFile) + " not valid.")

# record case info
caseinfo = {}
caseinfo['version'] = args.assayversion
caseinfo['name'] = args.name
caseinfo['mrn'] = args.mrn
caseinfo['DOB'] = args.DOB
caseinfo['accession'] = args.accession
caseinfo['specimen'] = args.specimen
caseinfo['casedir'] = args.dir
caseinfo['exception'] = args.exception
caseinfo['run_info_str'] = args.runinfostr
caseinfo['date'] = strftime("%Y-%m-%d_%H:%M:%S", gmtime())

if caseinfo['run_info_str'] == 'NONE' or caseinfo['run_info_str'] == 'null':
    caseinfo['runid'] = 'NONE'
    caseinfo['instrument'] = 'NONE'
    caseinfo['spec'] = 'NONE'
    caseinfo['flowcell'] = 'NONE'
else:
    run_info = caseinfo['run_info_str'].split(',')
    caseinfo['runid'] = run_info[0]
    caseinfo['instrument'] = ' '.join((run_info[1], 'Side', run_info[2]))
    caseinfo['spec'] = 'x'.join(run_info[-4:])
    caseinfo['flowcell'] = ' '.join((run_info[3], run_info[4], caseinfo['spec']))
    
if caseinfo['specimen'] == 'BM':
    caseinfo['specimen'] = 'Bone Marrow'
elif caseinfo['specimen'] == 'PB':
    caseinfo['specimen'] = 'Peripheral Blood'

#########################################
#
# Get thresholds and reference values for QC Metrics
#
#########################################

qcranges = {}
qcranges['qcdatafile'] = qcFile

with open(qcranges['qcdatafile'], 'r') as json_file:
    qcranges = json.load(json_file)

# update commandline thresholds/cutoffs if passed -- this overrides the ones in the QC JSON
for arg_name, arg_value in vars(args).items():
    if arg_name in qcranges.keys():
        qcranges[arg_name] = arg_value

# get nonsynonymous annotations
nonSynon = qcranges['PARAMETERS']['NONSYNONYMOUS_ANNOTATIONS']

# this is the minimum coverage for genes/exons
minCoverageLevel = 'pct_above_'+str(qcranges['PARAMETERS']['MinCoverage'])

dbcon = None
if args.variantdb is not None:
    qcranges['variantdb'] = args.variantdb
    dbcon = sqlite3.connect(qcranges['variantdb'])
else:
    qcranges['variantdb'] = None

#########################################
#
# Set up dataframes for all data
#
#########################################

# this is for dragen/umi metrics
qcdf = pd.DataFrame(columns=['metric','value','qcrange','flag','qcmetric'])

# this is for exon/gene coverage metrics
covDf = pd.DataFrame(columns=['Gene','Type','Region','Mean','Covered1','Covered2'])

# dataframe with small variants
variants = pd.DataFrame(columns=['type','filter','chrom','pos','ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations','coverage','altreads','vaf','dblookup'])

# dataframe with CNVs
svs = pd.DataFrame(columns=['category','type','chrom1','pos1','chrom2','pos2','length','csyntax','psyntax','genes','filters','id','abundance','Info','dblookup'])

# df for haplotect loci
haplotectlocidf = pd.DataFrame(columns=['chr','SNP1','SNP2','all11','all12','all21','all22','popn_counts','distance','total_count','sample_counts','contam_fraction'])


####
#
# Setup/open variant db
#

dbcon = None
if qcranges['variantdb'] is not None and 'POSITIVE' not in caseinfo['exception'].upper():
    # create database, if it doesnt exist
    if not os.path.exists(qcranges['variantdb']) or os.path.getsize(qcranges['variantdb']) == 0:
        print(f"Creating variant database {qcranges['variantdb']}")
        dbcon = sqlite3.connect(qcranges['variantdb'])
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),variants]).to_sql('variants', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),svs]).to_sql('svs', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),covDf]).to_sql('coverage', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),qcdf]).to_sql('qcdata', dbcon, if_exists='replace',index=False)
        pd.concat([pd.DataFrame(columns=caseinfo.keys()),haplotectlocidf]).to_sql('haplotectloci', dbcon, if_exists='replace',index=False)
    else:
        dbcon = sqlite3.connect(qcranges['variantdb'])


#########################################
#
# Get files from the case directory 
#
#########################################

vcffile = list(Path(caseinfo['casedir']).rglob('*.annotated.vcf.gz'))[0]
if not vcffile.is_file():
    sys.exit("VCF file " + str(vcffile) + " not valid.")

svvcffile = list(Path(caseinfo['casedir']).rglob('*.sv_annotated.vcf.gz'))[0]
if not svvcffile.is_file():
    sys.exit("SV VCF file " + str(svvcffile) + " not valid.")

cnvvcffile = list(Path(caseinfo['casedir']).rglob('*.cnv_annotated.vcf.gz'))[0]
if not cnvvcffile.is_file():
    sys.exit("CNV VCF file " + str(cnvvcffile) + " not valid.")

haplotect = list(Path(caseinfo['casedir']).rglob('*.haplotect.txt'))[0]
if not haplotect.is_file():
    sys.exit("Haplotect output " + str(haplotect) + " not valid.")

# dragen qc files
mappingMetrics = list(Path(caseinfo['casedir']).rglob('*.mapping_metrics.csv'))[0]
wgsMetrics = list(Path(caseinfo['casedir']).rglob('*.wgs_coverage_metrics.csv'))[0]
cnvMetrics = list(Path(caseinfo['casedir']).rglob('*.cnv_metrics.csv'))[0]
covReport = list(Path(caseinfo['casedir']).rglob('*.coverage_report.tsv'))[0]

if not mappingMetrics.is_file() or not wgsMetrics.is_file() or not cnvMetrics.is_file() or not covReport.is_file():
    sys.exit("DRAGEN metrics files not found.")

#########################################
#
# Collect QC metrics
#
#########################################

print("Collecting DRAGEN qc metrics...",file=sys.stderr)

# read in coverage metrics
df = pd.concat([pd.read_csv(mappingMetrics,sep=',',names=['group','readgroup','metric','value','percent']),
                pd.read_csv(wgsMetrics,sep=',',names=['group','readgroup','metric','value','percent']),
                pd.read_csv(cnvMetrics,sep=',',names=['group','readgroup','metric','value','percent'])],axis=0)

df['metric'] = df['group'] + ': ' + df['metric']
df = df.drop(columns='group')
dfpct = df[df['percent'].notna()].copy()
df = df.drop(columns='percent')
dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
dfpct['value'] = dfpct['percent']
dfpct = dfpct.drop(columns='percent')

# identify the ones we track
df['qcmetric'] = 0
df.loc[df['metric'].isin(qcranges['QC'].keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,df])
dfpct['qcmetric'] = 0
dfpct.loc[dfpct['metric'].isin(qcranges['QC'].keys()),'qcmetric'] = 1
qcdf = pd.concat([qcdf,dfpct])
pattern = r'\[\s*(\d+x): inf\)'
qcdf['metric'] = qcdf['metric'].str.replace(pattern, r'>\1', regex=True)
qcdf[['qcrange','flag']] = qcdf.apply(lambda r: checkQc(float(r['value']),qcranges['QC'][r['metric']]) if r['metric'] in qcranges['QC'].keys() else ['',''],axis=1,result_type='expand')
qcdf = qcdf[['metric','value','qcrange','flag','qcmetric']]

# update database with new results, if necessary
if dbcon is not None:
    dbUpdateVariants(dbcon,'qcdata',caseinfo,qcdf)

#########################################
#
# Collect coverage metrics
#
#########################################

with open(covReport, 'r') as f:
    lines = f.readlines()

# Remove '#' from the first line and use it as the header
header_line = lines[0].replace('#', '').strip()
covDf = pd.read_csv(covReport, header=None, skiprows=1, names=['Chromosome','Start','End','Gene','Info'] + header_line.split('\t')[5:],sep="\t")

# update database with new results, if necessary
if dbcon is not None:
    dbUpdateVariants(dbcon,'coverage',caseinfo,covDf)



geneCovDf = covDf[covDf['Info'].str.contains('CS_genes')]
svCovDf = covDf[covDf['Info'].str.contains('CS_svgenes')]

#
# get transcript info for genes and sv genes
#
ids = geneCovDf['Info'].str.split('|',expand=True).replace('\s\+\s\S+','',regex=True).loc[:,2:]
ids.columns = ['Transcript','Region','cdsStart','cdsEnd','strand']
geneTrx = pd.concat([geneCovDf,ids],axis=1)[['Gene','Transcript','cdsStart','cdsEnd','strand']].drop_duplicates().reset_index()
geneTrx['cdsStart'] = geneTrx['cdsStart'].astype(int)
geneTrx['cdsEnd'] = geneTrx['cdsEnd'].astype(int)

ids = svCovDf['Info'].str.split('|',expand=True).replace('\s\+\s\S+','',regex=True).loc[:,2:]
ids.columns = ['Transcript','Region','cdsStart','cdsEnd','strand']
svTrx = pd.concat([svCovDf,ids],axis=1)[['Gene','Transcript','cdsStart','cdsEnd','strand']].drop_duplicates()
svTrx['cdsStart'] = svTrx['cdsStart'].astype(int)
svTrx['cdsEnd'] = svTrx['cdsEnd'].astype(int)

# list of all genes and transcripts
knownTrx = pd.concat([svTrx[['Gene','Transcript']],geneTrx[['Gene','Transcript']]],axis=0,ignore_index=True).drop_duplicates()

covDf['Region'] = covDf['Info'].str.split('|').apply(lambda x: x[3] if len(x) > 2 else '')

geneList = covDf['Gene'].unique().tolist()

####################
#
# Haplotect
#
####################

haplotectlocidf = pd.read_csv(haplotect,sep='\t')
haplotectlocidf.columns = haplotectlocidf.columns.str.replace('#', '')
haplotectdf = pd.DataFrame([haplotectlocidf.iloc[-1,:-3].tolist()],columns=haplotectlocidf.iloc[-2,:-3].tolist())
haplotectdf.columns = haplotectdf.columns.str.replace('#', '')
haplotectlocidf = haplotectlocidf.iloc[:-2]

# update database with new results, if necessary
if dbcon is not None:
    dbUpdateVariants(dbcon,'haplotect',caseinfo,haplotectlocidf)


#########################################
#
# Get small variants
#
#########################################

print("Gathering gene variants...",file=sys.stderr)

vcf = VCF(vcffile)

# get VEP fields
vepFields = getVepFields(vcf)

# get variants
for variant in vcf:

    vartype = ''
    if len(variant.REF) == len(variant.ALT[0]):
        vartype = 'SNV'
    else:
        vartype = 'INDEL'

    varfilter = 'PASS'
    if variant.FILTER is not None:
        varfilter = variant.FILTER
 
    abundance = 'NA'
    totalReads = 'NA'
    variantReads = 'NA'

    abundance = round(variant.format('AF')[0][0] * 100,2)
    totalReads = variant.format("DP")[0][0]
    variantReads = variant.format("AD")[0][1]
        
    # get VEP annotation
    csq = variant.INFO['CSQ']
    
    if csq is None:
        sys.exit("No VEP fields")
    
    gene='NA'
    transcript='NA'
    csyntax='NA'
    psyntax='NA'
    consequence='NA'
    exon='NA'
    intron='NA'
    popmaf = 'NA'
    customannotation = 'NA'
    for i in variant.INFO['CSQ'].split(','):
        csq = i.split("|")

        # skip blacklisted variants
        if csq[vepFields['CHROMOSEQ_BLACKLIST']]:
            continue

        # get pop allele frequency. This is present for each transcript annotation, but is always the same 
        if csq[vepFields['MAX_AF']] != '':
            popmaf = float(csq[vepFields['MAX_AF']])

        # check if this is in the list of transcripts. only variants annotated with a known transcript will be reported
        if geneTrx['Transcript'].str.contains(csq[vepFields['Feature']]).any():
            transcript = csq[vepFields['Feature']]
            gene = csq[vepFields['SYMBOL']]
            consequence = csq[vepFields['Consequence']].split("&")[0]

            csyntax = csq[vepFields['HGVSc']].split(":")
            if len(csyntax) > 1:
                csyntax = csyntax[1]
            else:
                if csq[vepFields['STRAND']]!=geneTrx[geneTrx['Transcript']==transcript]['strand'].tolist()[0]:
                    sys.exit('strands dont match')

                if 'upstream' in consequence:
                    if csq[vepFields['STRAND']]==1:
                        # *  ---->                        
                        distance = geneTrx[geneTrx['Transcript']==transcript]['cdsStart'].min() - variant.POS
                        csyntax = "c.-"+str(distance)+variant.REF+'>'+csq[0]
                    else:
                        # <---- *
                        distance = variant.POS - geneTrx[geneTrx['Transcript']==transcript]['cdsEnd'].max()
                        csyntax = "c.-"+str(distance)+revcomp(variant.REF)+'>'+revcomp(csq[0])

                elif 'downstream' in consequence:
                    if csq[vepFields['STRAND']]==1:
                        # ---->  *                        
                        distance = variant.POS - geneTrx[geneTrx['Transcript']==transcript]['cdsEnd'].max()
                        csyntax = "c.+"+str(csq[vepFields['DISTANCE']])+variant.REF+'>'+csq[0]
                    else:
                        # *  <----
                        distance = geneTrx[geneTrx['Transcript']==transcript]['cdsEnd'].min() - variant.POS
                        csyntax = "c.+"+str(csq[vepFields['DISTANCE']])+revcomp(variant.REF)+'>'+revcomp(csq[0])
                else:
                    csyntax = 'noncoding'

            psyntax = csq[vepFields['HGVSp']].split(":")
            if len(psyntax) > 1:
                psyntax = convert_aa(psyntax[1])
                psyntax = re.sub("\%3D","=",psyntax)
            else:
                psyntax = csyntax
        
            impact = csq[vepFields['IMPACT']]
            exon = csq[vepFields['EXON']] or 'NA'
            intron = csq[vepFields['INTRON']] or 'NA'
            customannotation = csq[vepFields['Existing_variation']] or 'NA'    
            if 'CHROMOSEQ_DB' in vepFields.keys() and csq[vepFields['CHROMOSEQ_DB']]:
                customannotation = 'CSDB=' + str(csq[vepFields['CHROMOSEQ_DB']] or 0)        

    # convert pop maf to percent
    if popmaf!='NA':
        popmaf = round(float(popmaf)*100,3)

    # skip all variants > 0.1% and not in a hotspot 
    if popmaf!='NA' and popmaf >= float(args.popaf) and 'CHROMOSEQ' not in customannotation:
        continue

    # only include all variants <=0.1% and ns or specific noncoding variants 
    if consequence in nonSynon:
        variants = pd.concat([variants,pd.DataFrame([dict(zip(variants.columns,[vartype,varfilter,str(variant.CHROM),str(variant.POS),variant.REF,variant.ALT[0],gene,transcript,consequence,csyntax,psyntax,exon,str(popmaf) + '%',customannotation,str(totalReads),str(variantReads),str(abundance)+"%",'NA']))])])

# now query database to get counts of variants like this one.
if qcranges['variantdb'] is not None and variants.shape[0] > 0:
    queryPass = '''SELECT COUNT(*) FROM variants WHERE filter = "PASS" AND name != ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?'''
    queryFiltered = '''SELECT COUNT(*) FROM variants WHERE filter != "PASS" AND name != ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?'''
    variants['dblookup'] = variants.apply(lambda x: dbAnnotateVariants(dbcon,queryPass,queryFiltered,(caseinfo['name'],x['chrom'],x['pos'],x['ref'],x['alt'])),axis=1)

########################
#
# Get CNVs
#
########################

print("Gathering CNVs...",file=sys.stderr)

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

    category = 'CNV'
    if filter!='PASS':
        category = 'OtherSv'
    
    svs = pd.concat([svs,pd.DataFrame([dict(zip(svs.columns,[category,vartype,chr1,str(pos1),chr1,str(pos2),str(svlen),csyntax,psyntax,genestring,filter,str(variant.ID),str(abundance)+"%",'CN='+str(copynumber),'NA']))])])

########################
#
# Get SVs
# 
########################

print("Gathering SVs...",file=sys.stderr)

svvcf = VCF(svvcffile)

# get known fusion events
recurrentSvs = pd.read_csv(svGenes,sep=',', header=None)
recurrentSvs.columns = ['gene1','strand1','gene2','strand2']

passedvars = {}

for variant in svvcf:

    # save BNDs because we need both records for these
    if variant.INFO['SVTYPE']=='BND':
        passedvars[variant.ID] = variant

    else:

        filter = []

        vartype = variant.INFO.get('SVTYPE')

        # get gene 1 info
        chr = str(variant.CHROM)
        pos1 = str(variant.POS)
        pos2 = str(variant.INFO.get('END'))
        svlen = 'NA'
        if variant.INFO.get('SVLEN') is not None:
            svlen = str(abs(variant.INFO.get('SVLEN')))
        csyntax = ''
        psyntax = ''

        knownSvGenes = []
        if variant.INFO['KnownSvGenes'] is not None:
            knownSvGenes = knownSvGenes + variant.INFO['KnownSvGenes'].split(',')

        # get genes for this variant
        vepCsq = vepToTable(variant.INFO['CSQ'],svvcf.get_header_type('CSQ'))

        bands = vepCsq['cytobands'][0].split("&")
        bandstring = bands[0]
        if len(bands) > 1:
            bandstring = bands[0] + bands[-1]

        # indicate known genes by transcript id and then add exons if its a partial event.
        vepCsq['KnownTrx'] = 0
        vepCsq.loc[vepCsq['Feature'].isin(knownTrx['Transcript'].tolist()),'KnownTrx'] = 1
        vepCsq['KnownGene'] = 0
        vepCsq.loc[vepCsq['Feature'].isin(knownTrx['Gene'].tolist()),'KnownGene'] = 1

        vepCsq = vepCsq.sort_values(by=['KnownGene','KnownTrx','PICK'], ascending=[False,False,False]).drop_duplicates(subset='SYMBOL',keep='first')
        vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
        vepCsq['GeneImpact'] = vepCsq['Consequence'].apply(lambda r: int(len(set(r.split('&')) & set(nonSynon))>0))

        knownGeneDf = pd.concat([vepCsq[vepCsq['KnownGene']==1][['SYMBOL','GeneImpact','GeneEffect']],pd.DataFrame({'SYMBOL':knownSvGenes,'GeneEffect':''})]).sort_values(by=['SYMBOL','GeneImpact','GeneEffect'],ascending=[True,False,True],na_position='last').drop_duplicates(subset='SYMBOL',keep='first')
        # make a known gene string to report
        knownGeneString = 'None'
        if len(knownSvGenes)>0:
            knownGeneString = ','.join(knownGeneDf.apply(lambda x: x['SYMBOL'] + x['GeneEffect'],axis=1).tolist())

        # determine if this is a recurrent SV--if there are 1 or 2 known genes affected by this event and they match events in our list
        # since these are del/dup/ins events, the strands must be the same
        isRecurrentSv = False
        geneHits = knownGeneDf[knownGeneDf['GeneImpact']==1]['SYMBOL'].tolist()
        genePairHit = recurrentSvs[recurrentSvs.gene1.isin(geneHits) & recurrentSvs.gene2.isin(geneHits)]
        if genePairHit.shape[0] > 0 and genePairHit[genePairHit['strand1']==genePairHit['strand2']].shape[0] > 0:
            isRecurrentSv = True

        if vartype == 'DEL':
            csyntax = chr + ":g." + str(pos1) + "_" + str(pos2) + "del"
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                psyntax = "seq[GRCh38] -" + chr.replace('chr','')
                
            elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
                psyntax = "seq[GRCh38] -" + chr.replace('chr','')

            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "seq[GRCh38] del(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "seq[GRCh38] del(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
        
        elif vartype == 'DUP' or vartype == 'INS':
            csyntax = chr + ":g." + str(pos1) + "_" + str(pos2) + vartype.lower()
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                psyntax = "seq[GRCh38] +" + chr.replace('chr','')

            elif 'q11' in bands and 'qter' in bands and chr in ACROCENTRICS:
                psyntax = "seq[GRCh38] +" + chr.replace('chr','')

            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "seq[GRCh38] " + vartype.lower() + "(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "seq[GRCh38] " + vartype.lower() + "(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"

        abundance = 0.0
        PR = (0,0)
        SR = (0,0)            
        if variant.format("SR") is not None:
            SR = variant.format("SR")[0]
            
        if variant.format("PR")[0] is not None:                
            PR =  variant.format("PR")[0]

        abundance = (SR[1] + PR[1]) / (PR[0] + PR[1] + SR[0] + SR[1])

        infostring = 'PR_READS=' + str(PR[1]) + '/' + str(PR[0]+PR[1]) + ';SR_READS=' + str(SR[1]) + '/' + str(SR[0]+SR[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))

        if PR[1] < qcranges['PARAMETERS']['MinSvReads'] or SR[1] < qcranges['PARAMETERS']['MinSvReads']:
            filter = filter + ['MinSvReads']

        if abundance < qcranges['PARAMETERS']['MinSvAbundance']:
            filter = filter + ['MinSvAbundance']
            
        if variant.INFO.get('CONTIG') is None:
            filter = filter + ['Imprecise']

        if len(filter) == 0:
            filter = 'PASS'
        else:
            filter = ';'.join(filter)

        # categories and filter variants
        category = ''
        # if recurrent is true, then keep regardless
        if isRecurrentSv is True:
            category = 'Recurrent'

        # otherwise only keep PASS variants
        elif filter=='PASS':
            category = 'OtherSv'

        else:
            continue

        svs = pd.concat([svs,pd.DataFrame([dict(zip(svs.columns,[category,vartype,chr1,str(pos1),chr2,str(pos2),svlen,csyntax,psyntax,knownGeneString,filter,str(variant.ID),str(abundance)+"%",infostring,'NA']))])])

# now handle BNDs, which each have 2 entries.
# this includes translocations and inversions
alreadydone = set()

i  = 0
for v in passedvars.items():

    # get the first record
    variant = v[1]

    # skip if this is the mate and we already processed the record pair
    if variant.INFO.get('MATEID') in alreadydone or variant.INFO.get('MATEID') not in passedvars:
        continue

    # get the mate
    mate = passedvars[variant.INFO.get('MATEID')]

    vartype = variant.INFO.get('SVTYPE')

    # check to see if this is an inversion
    if variant.INFO.get('INV3') is not None or variant.INFO.get('INV5') is not None:
        vartype = 'INV'
    
    # get filter info
    filter = []
    if variant.FILTER is not None:
        filter = variant.FILTER.split(';')
    
    # get gene 1 info
    chr1 = str(variant.CHROM)
    pos1 = variant.POS
    genes1='NA'
    transcript1='NA'
    region1 = 'NA'
    strand1 = '+'
    bands1 = 'NA'

    chr2 = mate.CHROM
    pos2 = mate.POS     
    genes2='NA'
    transcript2='NA'
    region2 = 'NA'
    strand2 = '+'
    bands2 = 'NA'

    # get gene info for VARIANT
    # This gets each annota
    vepCsq = vepToTable(variant.INFO['CSQ'],svvcf.get_header_type('CSQ'))
    vepCsq['KnownTrx'] = 0
    vepCsq.loc[vepCsq['Feature'].isin(knownTrx['Transcript'].tolist()),'KnownTrx'] = 1
    vepCsq['KnownGene'] = 0
    vepCsq.loc[vepCsq['SYMBOL'].isin(knownTrx['Gene'].tolist()),'KnownGene'] = 1

    bands1 = vepCsq['cytobands'][0].split("&")[0]

    vepCsq = vepCsq.sort_values(by=['KnownGene','DISTANCE','KnownTrx','PICK'], ascending=[False,True,False,False]).drop_duplicates(subset='SYMBOL',keep='first')
    vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
    gene1Df= vepCsq[['SYMBOL','KnownGene','KnownTrx','DISTANCE','STRAND','GeneEffect']]

    if variant.INFO.get('KnownSvGenes') is not None:
        gene1Df = pd.concat([gene1Df,pd.merge(pd.DataFrame({'SYMBOL':variant.INFO['KnownSvGenes'].split(','),'KnownGene':1}),
                           pd.concat([recurrentSvs[['gene1','strand1']].rename(columns={'gene1':'SYMBOL','strand1':'STRAND'}),recurrentSvs[['gene2','strand2']].rename(columns={'gene2':'SYMBOL','strand2':'STRAND'})],axis=0).drop_duplicates(),on='SYMBOL')],axis=0)

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

# now query database to get counts of variants like this one.
if qcranges['variantdb'] is not None and svs.shape[0] > 0:
    queryPass = '''SELECT COUNT(*) FROM svs WHERE filter = "PASS" AND name != ? AND gene1 = ? AND region1 = ? AND gene2 = ? AND region2 = ?'''
    queryFiltered = '''SELECT COUNT(*) FROM svs WHERE filter != "PASS" AND name != ? AND gene1 = ? AND region1 = ? AND gene2 = ? AND region2 = ?'''
    svs['dblookup'] = svs.apply(lambda x: dbAnnotateVariants(dbcon,queryPass,queryFiltered,(caseinfo['name'],x['gene1'],x['region1'],x['gene2'],x['region2'])),axis=1)

# close database
if qcranges['variantdb'] is not None:
    dbcon.close()

########################
#
# Start report
#
########################

print("Starting report...",file=sys.stderr)

# make dict for report and redirect output for text report
jsonout = {'CASEINFO':{},'VARIANTS':{},'QC':{}}

outfile = open(caseinfo['name'] + ".chromoseq.txt", "w")
sys.stdout = outfile

print("ChromoSeq Report for " + caseinfo['name'] + " ---- Generated on: " + caseinfo['date'] + "\n")

print("*** CHROMOSEQ CASE INFORMATION ***\n")
print("MRN:\t" + caseinfo['mrn'])
print("ACCESSION:\t" + caseinfo['accession'])
print("SPECIMEN TYPE:\t" + caseinfo['specimen'])
print("DOB:\t" + caseinfo['DOB'])
print("RUNID:\t" + caseinfo['runid'])
print("INSTRUMENT:\t" + caseinfo['instrument'])
print("FLOWCELL:\t" + caseinfo['flowcell'])
if (caseinfo['exception'] != 'NONE'):
    caseinfo['exception'] = caseinfo['exception'] + "\t(!)"
print("EXCEPTIONS:\t" + caseinfo['exception'])
jsonout['CASEINFO'] = caseinfo

print()

print("*** SEQUENCING QC ***\n")

print(qcdf[qcdf['qcmetric']==1].drop(columns='qcmetric').to_csv(sep='\t',header=False, index=False))
jsonout['QC']['ASSAY'] = qcdf[qcdf['qcmetric']==1].drop(columns='qcmetric').to_dict('split')
del jsonout['QC']['ASSAY']['index']

print()

print("*** PLOIDY AND PURITY ***\n")

print('Purity:\t' + str(qcdf[qcdf['metric'].str.contains('tumor purity')]['value'].tolist()[0]))

ploidy = float(qcdf[qcdf['metric'].str.contains('Overall ploidy')]['value'].tolist()[0])
chromosomes = str(int(round(ploidy * 23,0)))
print('Ploidy:\t' + chromosomes + ' (' + str(ploidy) + ')')

print()

print("*** COPY NUMBER ALTERATIONS ***\n")

if svs[svs['category']=='CNV'].shape[0] > 0:
    print(svs[svs['category']=='CNV'].drop(columns='category').to_csv(sep='\t',header=True, index=False))
else:
    print("None Detected\n")

jsonout['VARIANTS']['CNV'] = svs[svs['category']=='CNV'].drop(columns='category').to_dict('split')
del jsonout['VARIANTS']['CNV']['index']

print("*** RECURRENT STRUCTURAL VARIANTS ***\n")

if svs[svs['category']=='Recurrent'].shape[0] > 0:
    print(svs[svs['category']=='Recurrent'].drop(columns='category').to_csv(sep='\t',header=True, index=False))    
else:
    print("None Detected\n")

jsonout['VARIANTS']['RECURRENTSV'] = svs[svs['category']=='Recurrent'].drop(columns='category').to_dict('split')
del jsonout['VARIANTS']['RECURRENTSV']['index']

print("*** GENE MUTATIONS ***\n")

if variants[variants['filter']=='PASS'].shape[0] > 0:
    print(variants[variants['filter']=='PASS'].to_csv(sep='\t',header=True, index=False))
else:
    print("None Detected\n")

jsonout['VARIANTS']['PASS'] = variants[variants['filter']=='PASS'].to_dict('split')
del jsonout['VARIANTS']['PASS']['index']

print("*** FILTERED GENE MUTATIONS ***\n")

if variants[variants['filter']!='PASS'].shape[0] > 0:
    print(variants[variants['filter']!='PASS'].to_csv(sep='\t',header=True, index=False))
else:
    print("None Detected\n")

jsonout['VARIANTS']['FILTERED'] = variants[variants['filter']!='PASS'].to_dict('split')
del jsonout['VARIANTS']['FILTERED']['index']

print("*** NOVEL/FILTERED SV AND CNV VARIANTS ***\n")

if svs[svs['category']=='OtherSv'].shape[0] > 0:
    print(svs[svs['category']=='OtherSv'].drop(columns='category').to_csv(sep='\t',header=True, index=False))    
else:
    print("None Detected\n")

jsonout['VARIANTS']['OTHERSV'] = svs[svs['category']=='OtherSv'].drop(columns='category').to_dict('split')
del jsonout['VARIANTS']['OTHERSV']['index']

print()

print("*** FAILED EXONS ***\n")

xdf = covDf[(covDf['Info'].str.contains('CS_genes'))].copy()
xdf['QC'] = ''
xdf.loc[xdf[minCoverageLevel]<qcranges['PARAMETERS']['MinFractionCoverage'],'QC'] = '(!)'
if xdf[xdf['QC']!=''].shape[0] > 0:
    print(xdf[xdf['QC']!=''][['Gene','Region','mean_cvg',minCoverageLevel,'QC']].to_csv(sep='\t',header=True, index=False,float_format='%.1f'))
else:
    print("No failed exons\n")

# only report failed exons
jsonout['QC']['FAILED EXON QC'] = xdf.to_dict('split')
del jsonout['QC']['FAILED EXON QC']['index']

print("*** FAILED GENE QC ***\n")

xdf = covDf[(covDf['Info'].str.contains('CS_svgenes'))].copy()
xdf['QC'] = ''
xdf.loc[xdf[minCoverageLevel]<qcranges['PARAMETERS']['MinFractionCoverage'],'QC'] = '(!)'
if xdf[xdf['QC']!=''].shape[0] > 0:
    print(xdf[xdf['QC']!=''][['Gene','Region','mean_cvg',minCoverageLevel,'QC']].to_csv(sep='\t',header=True, index=False,float_format='%.1f'))
else:
    print("No failed genes\n")

jsonout['QC']['GENE QC'] = xdf.to_dict('split')
del jsonout['QC']['GENE QC']['index']

print("*** Haplotect Contamination Estimate ***\n")

print(haplotectdf.iloc[:,1:].to_csv(sep='\t',header=True, index=False))
print(haplotectlocidf.to_csv(sep='\t',header=True, index=False))

jsonout['QC']['HAPLOTECT SUMMARY'] = haplotectdf.iloc[:,1:].to_dict('split')
del jsonout['QC']['HAPLOTECT SUMMARY']['index']

jsonout['QC']['HAPLOTECT LOCI'] = haplotectlocidf.iloc[:,1:].to_dict('split')
del jsonout['QC']['HAPLOTECT LOCI']['index']

print("*** ChromoSeq Assay Version " + str(caseinfo["version"]) + ") ***\n")

print(qcranges['PARAMETERS']["DISCLAIMER"])

jsonout['QC']['PARAMETERS'] = qcranges['PARAMETERS']
jsonout['QC']['QCINFO'] = qcranges['QC']

original_stdout = sys.stdout
outfile.close()

# dump json
j = open(caseinfo['name'] + ".chromoseq.json", "w")
json.dump(jsonout,j,indent=" ")
j.close()
