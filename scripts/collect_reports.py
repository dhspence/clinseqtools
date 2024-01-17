import argparse
import pandas as pd
import sys, os

parser = argparse.ArgumentParser(description="Collect ChromoSeq Reports.")
parser.add_argument("files", type=str, nargs='+', help="List of input CSV files.")
parser.add_argument("-o","--outfile", type=str, help="Output file. If not provided, prints to stdout.")
parser.add_argument("-t","--tables", type=str,default=None, help="Print individual tables for variants, svs, and qc.")

args = parser.parse_args()

variants = pd.DataFrame()
svs = pd.DataFrame()
qc = pd.DataFrame()
report = pd.DataFrame()

for file in args.files:
    if not os.path.exists(file):
        print(f"Error: File '{file}' does not exist.")
        continue

    rep = pd.read_json(file)

    caseInfo = pd.DataFrame(data=rep['CASEINFO'].dropna()).to_dict()['CASEINFO']

    variants = pd.concat([variants,pd.concat([pd.DataFrame(data=rep['VARIANTS'].dropna().to_dict()['PASS']['data'],
                                   columns=rep['VARIANTS'].dropna().to_dict()['PASS']['columns']).assign(category='VARIANTS').assign(name=caseInfo['name']),
                      pd.DataFrame(data=rep['VARIANTS'].dropna().to_dict()['FILTERED']['data'],
                                   columns=rep['VARIANTS'].dropna().to_dict()['FILTERED']['columns']).assign(category='FILTERED').assign(name=caseInfo['name'])])],ignore_index=True,axis=0)
    svs = pd.concat([svs,pd.concat([pd.DataFrame(data=rep['VARIANTS'].dropna().to_dict()['CNV']['data'],
                    columns=rep['VARIANTS'].dropna().to_dict()['CNV']['columns']).assign(category='CNV').assign(name=caseInfo['name']),
                 pd.DataFrame(data=rep['VARIANTS'].dropna().to_dict()['RECURRENTSV']['data'],
                   columns=rep['VARIANTS'].dropna().to_dict()['RECURRENTSV']['columns']).assign(category='RECURRENTSV').assign(name=caseInfo['name']),
                 pd.DataFrame(data=rep['VARIANTS'].dropna().to_dict()['OTHERSV']['data'],
                              columns=rep['VARIANTS'].dropna().to_dict()['OTHERSV']['columns']).assign(category='OTHERSV').assign(name=caseInfo['name'])])],ignore_index=True,axis=0)

    qc = pd.concat([qc,pd.concat([pd.DataFrame(data=rep['CASEINFO'].dropna()).transpose().reset_index().drop('index',axis=1),pd.DataFrame(data=rep['QC']['ASSAY']['data'],columns=rep['QC']['ASSAY']['columns'])[['metric','value']].set_index('metric').transpose().reset_index().drop('index',axis=1)],axis=1)],ignore_index=True,axis=0)


outdf = pd.merge(qc,variants[variants['category']=='VARIANTS'].assign(genevariants=variants['gene']+":"+variants['psyntax']+"["+variants['vaf']+"]").groupby('name')['genevariants'].agg(lambda x: ','.join(x)).reset_index(),on='name',how='left')
outdf = pd.merge(outdf,svs[svs['category']=='RECURRENTSV'].groupby('name')['genes'].agg(lambda x: ','.join(x)).reset_index().rename(columns={'genes':'recurrentsv'}),on='name',how='left')
outdf = pd.merge(outdf,svs[svs['category']=='CNV'].assign(cnvs=svs[svs['category']=='CNV']['psyntax'].str.replace('seq[GRCh38] ','')+"["+svs[svs['category']=='CNV']['genes']+"]").groupby('name')['cnvs'].agg(lambda x: ','.join(x)).reset_index(),on='name',how='left')
outdf = pd.merge(outdf,svs[svs['category']=='OTHERSV'].assign(othersvs=svs[svs['category']=='OTHERSV']['psyntax'].str.replace('seq[GRCh38] ','')+"["+svs[svs['category']=='OTHERSV']['genes']+"]").groupby('name')['othersvs'].agg(lambda x: ','.join(x)).reset_index(),on='name',how='left')

if args.tables is not None:
    variants.to_csv(args.tables+'_variants.tsv', sep="\t", index=False)
    svs.to_csv(args.tables+'_svs.tsv', sep="\t", index=False)
    qc.to_csv(args.tables+'_qc.tsv', sep="\t", index=False)

if args.outfile:
    outdf.to_csv(args.outfile, sep='\t', index=False)
else:
    outdf.to_csv(sys.stdout, sep='\t', index=False)

