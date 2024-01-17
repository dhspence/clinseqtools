#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import json, sys, os, csv, argparse, gzip

__version__ = '1.0.0'

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

def make_runinfo(read1path,read2path):

    # Open the gzipped file and read the first two lines
    readName = ''
    seq1 = ''
    seq2 = ''
    with gzip.open(read1path, 'rt') as gz_file:
        readName = gz_file.readline().strip()  # First line
        seq1 = gz_file.readline().strip()       # Second line

    # Open the gzipped file and read the first two lines
    with gzip.open(read2path, 'rt') as gz_file:
        readName = gz_file.readline().strip()  # First line
        seq2 = gz_file.readline().strip()       # Second line

    parts = readName.split(':')
    index1 = 'UNKNOWN'
    index2 = 'UNKNOWN'
    
    index1, index2 = readName.split(' ').split(':').split('+')
    # Extract required parts and format them
    runinfo = {'RunID':f"RUN_{parts[0][1:]}_{str(int(parts[1])).zfill(4)}_{parts[2]}",
                'Flowcell':parts[2],
                'Lane':parts[3],
                'Instrument':parts[0][1:],
                'Read1Cycles':len(seq1),
                "Index1Cycles": 'UNKNOWN',
                "Index1Reverse": "UNKNOWN",
                "Index2Cycles": "UNKNOWN",
                "Index2Reverse": "UNKNOWN",
                "Read2Cycles": len(seq2),
                "Index1":index1,
                "Index2":index2
    }

    return runinfo


def main():

    parser = argparse.ArgumentParser(description="Prepare fastq_list file from a Dragen demux directory")
    parser.add_argument("sample_id", type=str, help="Sample ID")
    parser.add_argument("read1", type=checkfile, help="Path to demux directory")
    parser.add_argument("read2", type=checkfile, help="Path to demux directory")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    id = args.sample_id
    read1 = args.read1
    read2 = args.read2

    info = make_runinfo(read1,read2)

    if info and 'Flowcell' in info:

        fastqlist_file = 'fastq_list.csv'
        with open(fastqlist_file, 'w', newline='') as cfile:
            writer = csv.writer(cfile)
            writer.writerow('RGID,RGSM,RGLB,Lane,Read1File,Read2File')
            rgid = '.'.join([info['Flowcell'],info['Index1'],info['Index2'],info['Lane']])
            rglb = '.'.join([id,info['Index1'],info['Index2']])
            writer.writerow([rgid,id,rglb,info['Lane'],os.path.basename(read1),os.path.basename(read2)])

        del info['Index1']
        del info['Index2']
        del info['Lane']
        csv_file_name = info['Flowcell'] + '.run_info.csv'
        with open(csv_file_name, 'w', newline='') as cfile:
            writer = csv.writer(cfile)
            writer.writerow(info.keys())
            writer.writerow(info.values())

    else:
        print("Error: Could not retrieve information or missing Flowcell ID.")



if __name__ == "__main__":
    main()

