#!/usr/bin/env python3

import pysam
import pandas as pd
import subprocess, sys, os
import numpy as np
import argparse
import gzip

__version__ = '1.0.0'

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def fileexists(file_path):
    """Check if a file exists at the given path."""
    if os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The ouytfile {file_path} exists!")
    return file_path

def main():
    parser = argparse.ArgumentParser(description="Coverage analysis using samtools depth.")
    parser.add_argument("bed_file", type=checkfile, help="Path to the bed file (plain text or gzipped).")
    parser.add_argument("cram_file", type=checkfile, help="Path to the cram/bam file.")
    parser.add_argument("-o","--out", type=fileexists, help="Out file")
    parser.add_argument("-r","--reference", type=checkfile, help="Reference fasta")    
    parser.add_argument("-q","--minmapqual", type=int, default=1, help="Minimum mapping quality")
    parser.add_argument("-Q","--minbasequal", type=int, default=13, help="Minimum base quality")
    parser.add_argument("-c","--coverage_values", default="10,20,40,60,100,250,500,1000,1250,2500,4000", 
                        help="Comma-separated coverage values (default: 10,20,40,60,100,250,500,1000,1250,2500,4000).")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()
    
    cov_values = [int(x) for x in args.coverage_values.split(",")]
    minMapQual = args.minmapqual
    minBaseQual = args.minbasequal

    if args.out:
        outfile = open(args.out, "w")
    else:
        outfile = sys.stdout
    
    header = ["#chrom", "start", "end", "gene", "info", "total_cvg", "mean_cvg", "Q1_cvg", "median_cvg", "Q3_cvg", "min_cvg", "max_cvg"]
    header += [f"pct_above_{val}" for val in cov_values]
    print("\t".join(header),file=outfile)

    # Extract regions from bed file
    with (gzip.open(args.bed_file, 'rt') if args.bed_file.endswith('.gz') else open(args.bed_file, 'r')) as bed:
        for line in bed:
            chrom, start, end, gene, info = line.strip().split()[:5]

            result  = pysam.depth("--reference",args.reference,"-q",str(minMapQual),"-Q",str(minBaseQual),"-s", "-r", f"{chrom}:{start}-{end}",args.cram_file,catch_stdout=True)
            
            data = [line.split() for line in result.splitlines()]
            df = pd.DataFrame(data, columns=['chrom', 'position', 'coverage'])
            df['coverage'] = df['coverage'].astype(int)

            sum_cov = df['coverage'].sum()
            mean_cov = round(df['coverage'].mean(), 1)
            q1 = round(df['coverage'].quantile(0.25), 1)
            median_cov = round(df['coverage'].median(), 1)
            q3 = round(df['coverage'].quantile(0.75), 1)
            min_cov = df['coverage'].min()
            max_cov = df['coverage'].max()
            perc_above = [round((df['coverage'] > val).mean() * 100, 1) for val in cov_values]

            print(chrom, start, end, gene, info, sum_cov, mean_cov, q1, median_cov, q3, min_cov, max_cov, *perc_above, sep="\t",file=outfile)

    outfile.close()

if __name__ == "__main__":
    main()
