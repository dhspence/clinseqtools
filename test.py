#!/usr/bin/env python3

import sys
import pandas as pd
from clinseqtools import ClinSeqTools


def main():

    myFile = 'example_data/TWDY-ChromoSeq-Validation-M1525590-lib1.annotated.vcf.gz' #sys.argv[1]
    myTranscripts = 'example_data/TWDY-ChromoSeq-Validation-M1525590-lib1.coverage_report.tsv' #sys.argv[1]

    myCase = ClinSeqTools()
    #myCase.parse_coverage_report(file=myTranscripts)
    myCase.parse_variant_file(file=myFile)
    print("made it")

if __name__ == '__main__':
    main()
    
        


