#!/usr/bin/env python3

import sys
import pandas as pd
from clinseqtools import ClinSeqTools


def main():

    myFile = 'TWDY-ChromoSeq-Validation-M1525590-lib1.annotated.vcf.gz' #sys.argv[1]

    myTranscripts = 'TWDY-ChromoSeq-Validation-M1525590-lib1.coverage_report.tsv' #sys.argv[1]
    with open(myTranscripts, 'r') as f:
        lines = f.readlines()

    # Remove '#' from the first line and use it as the header
    header_line = lines[0].replace('#', '').strip()
    covDf = pd.read_csv(myTranscripts, header=None, skiprows=1, names=['chromosome','start','end','gene','info'] + header_line.split('\t')[5:],sep="\t")
    geneCovDf = covDf[covDf['info'].str.contains('CS_genes')]
    ids = geneCovDf['info'].str.split('|',expand=True).replace('\s\+\s\S+','',regex=True).loc[:,2:]
    ids.columns = ['transcript','region','cdsstart','cdsend','strand']
    geneTrx = pd.concat([geneCovDf,ids],axis=1)[['gene','transcript','cdsstart','cdsend','strand']].drop_duplicates()
    geneTrx[['cdsstart','cdsend','strand']] = geneTrx[['cdsstart','cdsend','strand']].astype(int)

    myCase = ClinSeqTools()
    myCase.parse_variant_file(file=myFile,transcripts=geneTrx)
    print("made it")

if __name__ == '__main__':
    main()
    
        


