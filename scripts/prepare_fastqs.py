#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import json, sys, os, csv, argparse, gzip
from collections import Counter
import pandas as pd

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

def make_runinfo_from_reads(read1path,read2path):

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
    indexes = readName.split(' ')
    if len(indexes) > 1:
        indexlist = []
        with gzip.open(read1path, 'rt') as file:
            for i, line in enumerate(file):
                if i % 4 == 0:  # Read names are on every 4th line starting from 0
                    read_name = line.strip()  # Remove leading and trailing whitespaces
                    indexes = read_name.split(' ')
                    if len(indexes) > 1:
                        index = indexes[1].split(':')[-1]
                        indexlist = indexlist + [index]
        
                if i == 399:  # Stop after reading the first 100 entries (0 to 399)
                    break

        counter = Counter(indexlist)
        most_common_index, _ = counter.most_common(1)[0]
        index1, index2 = most_common_index.split('+')

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

def make_runinfo_from_demuxpath(file_path):

    runinfo_path = os.path.join(file_path, "Reports", "RunInfo.xml")

    if not os.path.exists(runinfo_path):
        print(f"Error: '{runinfo_path}' does not exist in {file_path}.")
        sys.exit(1)

    # Parse the XML file
    tree = ET.parse(runinfo_path)
    root = tree.getroot()

    # Initialize dictionary to store the required data
    run_info = {}

    # Extract the Run ID, Flowcell, and Instrument
    run_dat = root.find('Run')
    if run_dat is not None:
        run_info['RunID'] = run_dat.get('Id')
        run_info['Flowcell'] = run_dat.find('Flowcell').text
        run_info['Instrument'] = run_dat.find('Instrument').text

    # Extract NumCycles for each read
    reads = run_dat.find('Reads')
    for i, read in enumerate(reads.findall('Read'), start=1):
        num_cycles = read.get('NumCycles')
        if num_cycles is not None:
            if int(i) == 1:
                run_info[f'Read1Cycles'] = num_cycles
            elif int(i) == 2:
                run_info[f'Index1Cycles'] = num_cycles
                run_info[f'Index1Reverse'] = read.get("IsReverseComplement")
            elif int(i) == 3:
                run_info[f'Index2Cycles'] = num_cycles
                run_info[f'Index2Reverse'] = read.get("IsReverseComplement")
            elif int(i) == 4:
                run_info[f'Read2Cycles'] = num_cycles

    runpars_path = os.path.join(file_path, "Reports", "RunParameters.xml")

    if os.path.exists(runpars_path):
        tree = ET.parse(runpars_path)
        root = tree.getroot()

        run_info['FlowCellType'] = root.find('FlowCellType').text
        run_info['InstrumentType'] = root.find('InstrumentType').text
        run_info['Side'] = root.find('Side').text
        # fine the serial number of the flowcell:
        for el in root.findall('ConsumableInfo')[0].findall('ConsumableInfo'):
            el2 = el.find('Type')
            if el2 is not None:
                if el2.text == 'FlowCell':
                    run_info['FlowCellType'] = el.find('Mode').text
                    run_info['FlowCellLotNumber'] = el.find('LotNumber').text
                elif el2.text == 'Reagent':
                    run_info['ReagentLotNumber'] = el.find('LotNumber').text

    return run_info


def main():

    parser = argparse.ArgumentParser(description="Prepare fastq_list file from a Dragen demux directory")
    parser.add_argument('-i', '--id', type=str, required=True, help="Sample ID")
    parser.add_argument('-1','--read1', type=checkfile, help="Path to read1")
    parser.add_argument('-2','--read2', type=checkfile, help="Path to read2")
    parser.add_argument('-f','--fastqlist', type=checkfile, help="Path to fastq_list.csv")
    parser.add_argument('-d','--demuxpath', type=checkfile, help="Path to Dragen demux directory")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    id = args.id

    # if read1 and read2 are supplied
    if args.read1 and args.read2:
        read1 = os.path.realpath(args.read1)
        read2 = os.path.realpath(args.read2)
        runinfo = make_runinfo_from_reads(read1,read2)
        
        fastqlist_file = 'fastq_list.csv'
        with open(fastqlist_file, 'w', newline='') as cfile:
            writer = csv.writer(cfile)
            writer.writerow(['RGID','RGSM','RGLB','Lane','Read1File','Read2File'])
            rgid = '.'.join([runinfo['Flowcell'],runinfo['Index1'],runinfo['Index2'],runinfo['Lane']])
            rglb = '.'.join([id,runinfo['Index1'],runinfo['Index2']])
            writer.writerow([rgid,id,rglb,runinfo['Lane'],read1,read2])

        del runinfo['Index1']
        del runinfo['Index2']
        del runinfo['Lane']
        csv_file_name = id + '.run_info.csv'
        with open(csv_file_name, 'w', newline='') as cfile:
            writer = csv.writer(cfile)
            writer.writerow(runinfo.keys())
            writer.writerow(runinfo.values())

    # if a fastq list is supplied
    elif args.fastqlist or args.demuxpath:
        runinfoDf = pd.DataFrame()
        fqlist = pd.DataFrame()
        if args.fastqlist:
            fqlist = pd.read_csv(args.fastqlist)
            fqlist = fqlist[fqlist['RGSM']==id]
            runinfo = []
            # Iterate through each row in the DataFrame
            for index, row in fqlist.iterrows():
                # Extract values from the specified columns
                read1 = row['Read1File']
                read2 = row['Read2File']
                # Apply the make_info function and add the result to the list
                runinfo = runinfo + [make_runinfo_from_reads(row['Read1File'], row['Read2File'])]

            # Convert the list of dictionaries to a DataFrame
            runinfoDf = pd.DataFrame(runinfo).drop_duplicates()

        elif args.demuxpath:
            runinfo = make_runinfo_from_demuxpath(args.demuxpath)
            fastqlistFile = os.path.join(args.demuxpath,'Reports','fastq_list.csv')
            if not os.path.exists(fastqlistFile):
                print(f"Fastq list not found in demux path: {args.demuxpath}")
                sys.exit(1) 

            fqlist = pd.read_csv(fastqlistFile)
            fqlist = fqlist[fqlist['RGSM']==id]

            runinfoDf = pd.DataFrame([make_runinfo_from_demuxpath(args.demuxpath)])

        fqlist.to_csv('fastq_list.csv',index=False)
        runinfoDf.to_csv(id + '.run_info.csv',index=False)

if __name__ == "__main__":
    main()

