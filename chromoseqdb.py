import os, sys, re
import shutil
import argparse
import pandas as pd
from datetime import datetime

from sqlalchemy import (create_engine, inspect, MetaData, Column, Integer, String, Float,
                        DateTime, Date)
from sqlalchemy.orm import declarative_base, sessionmaker

Base = declarative_base()

# Define ORM classes
class cases(Base):
    __tablename__ = 'cases'
    key = Column(Integer, primary_key=True, autoincrement=True)
    caseid = Column(String)
    date = Column(String)
    version = Column(String)
    mrn = Column(Integer)
    accession = Column(String)
    specimen_type = Column(String)
    dob = Column(String)
    sex = Column(String)
    runid = Column(String)
    instrument = Column(String)
    flowcell = Column(String)

class qc(Base):
    __tablename__ = 'qc'
    key = Column(Integer, primary_key=True, autoincrement=True)
    caseid = Column(String)
    date = Column(String)
    version = Column(String)
    total_reads = Column(Integer)
    total_bases = Column(Integer)
    insert_size = Column(Float)
    coverage = Column(Float)
    duplicate_rate = Column(Float)
    mapped_reads = Column(Float)
    mismatched_bases_R1 = Column(Float)
    mismatched_bases_R2 = Column(Float)
    q30_bases_R1 = Column(Float)
    q30_bases_R2 = Column(Float)
    haplotect_sites = Column(Integer)
    haplotect_informative_sites = Column(Integer)
    haplotect_coverage = Column(Float)
    haplotect_total_reads = Column(Integer)
    haplotect_contaminating_reads = Column(Integer)
    haplotect_contamination_fraction = Column(Float)
    haplotect_contamination_estimate = Column(Float)
    haplotect_contamination_confidence = Column(String)
    purity = Column(String)
    ploidy = Column(String)

class svs(Base):
    __tablename__ = 'svs'
    key = Column(Integer, primary_key=True, autoincrement=True)
    caseid = Column(String)
    date = Column(String)
    version = Column(String)
    category = Column(String)
    type = Column(String)
    chrom1 = Column(String)
    pos1 = Column(Integer)
    chrom2 = Column(String)
    pos2 = Column(Integer)
    length = Column(Integer)
    csyntax = Column(String) 
    psyntax = Column(String) 
    genes = Column(String)   
    filters = Column(String) 
    id = Column(String)
    abundance = Column(String)
    info = Column(String)


class variants(Base):
    __tablename__ = 'variants'
    key = Column(Integer, primary_key=True, autoincrement=True)
    caseid = Column(String)
    date = Column(String)
    version = Column(String)
    category = Column(String)
    type = Column(String)
    filters = Column(String)
    chrom = Column(String)
    pos = Column(Integer)
    ref = Column(String)
    alt = Column(String)
    gene = Column(String)
    transcript = Column(String)
    consequence = Column(String)
    csyntax = Column(String)
    psyntax = Column(String)
    exon = Column(String)
    popaf = Column(Float)
    annotations = Column(String)
    coverage = Column(Integer)
    altreads = Column(Integer)
    vaf = Column(Float)

class history(Base):
    __tablename__ = 'history'
    key = Column(Integer, primary_key=True, autoincrement=True)
    date = Column(DateTime)
    records = Column(Integer)


def create_db(filepath):
    print(f'Creating ChromoSeq database {filepath}',file=sys.stderr)

    engine = create_engine(f'sqlite:///{filepath}')
    Base.metadata.create_all(engine)

    # create a history record
    Session = sessionmaker(bind=engine)
    session = Session()
    new_history_entry = history(date=datetime.utcnow(), records=0)
    session.add(new_history_entry)
    session.commit()
    session.close()

def add_records_from_dataframe(df, table_class, engine):
    # Create a session
    Session = sessionmaker(bind=engine)
    session = Session()

    # 1. Get the DataFrame column names
    df_columns = set(df.columns.tolist())

    # 2. Get the database table column names
    inspector = inspect(engine)
    table_columns = inspector.get_columns(table_class.__tablename__)
    db_columns = set([column['name'] for column in table_columns])

    # 3. Compare columns while ignoring the primary key column 'key'
    db_columns.discard('key')

    if df_columns != db_columns:
        session.close()
        raise ValueError("DataFrame columns and database table columns do not match!")

    # delete any records with the same caseid. We assume this is a repeat update
    unique_caseids = df['caseid'].unique()
    for caseid in unique_caseids:
        session.query(table_class).filter(table_class.caseid == caseid).delete()

    session.commit()

    df.to_sql(table_class.__tablename__, con=engine, if_exists='append', index=False)
    session.close()

def query_table(table_class, engine):
    Session = sessionmaker(bind=engine)
    session = Session()
    
    records = session.query(table_class).all()
    column_names = [column.name for column in table_class.__table__.columns]
    primary_keys = [column.name for column in table_class.__table__.columns if column.primary_key]

    df = pd.DataFrame([[getattr(record, column) for column in column_names] for record in records], columns=column_names)
    df = df.drop(primary_keys, axis=1)

    session.close()
    
    return df

def backup_if_needed(db, engine, forcebackup):
    Session = sessionmaker(bind=engine)
    session = Session()
    last_date = session.query(history.date).order_by(history.date.desc()).first()

    # if no record or older than 30 days or force is true then add a record
    if (datetime.now() - last_date[0]).days > 30 or forcebackup:
        print(f'Backup up database {db}',file=sys.stderr)
       
        records = session.query(cases).all()
        record_count = len(records)
        new_history_entry = history(date=datetime.now(), records=record_count)
        session.add(new_history_entry)
        session.commit()
        baseName = os.path.splitext(os.path.basename(db))[0]
        backup_name = os.path.join(os.path.dirname(db),f"{baseName}_{datetime.now().strftime('%Y%m%d')}.db")
        shutil.copy(db, backup_name)

    session.close()

def parse_cs_report_text(filename):

    categories = {'COPY NUMBER ALTERATIONS':'CNV', 
                  'RECURRENT TRANSLOCATIONS':'RECURRENTSV',
                  'OTHER STRUCTURAL VARIANTS':'OTHERSV',
                  'NOVEL/FILTERED SV AND CNV VARIANTS':'OTHERSV',
                  'GENE MUTATIONS':'PASS',
                  'FILTERED GENE MUTATIONS':'FILTERED'}
    
    with open(filename, 'r') as f:
        content = f.read()

    # Extract report id and date
    match = re.match(r'ChromoSeq Report for (.+) ---- Generated on: (\d{4}-\d{2}-\d{2}.\d{2}:\d{2}:\d{2})', content)
    if match:
        report_id, date = match.groups()
    else:
        print("Header format is incorrect.")
        return

    # Split the content into sections based on ***
    sections = re.split(r'\*\*\* ([A-Za-z\s]+) \*\*\*', content)[1:]

    csVersion = 'beta'
    if 'Version' in sections[-1]:
        csVersion = re.split(r'\*\*\*',sections[-1])
        if len(csVersion) < 3:
            csVersion = 'beta'
        else:    
            csVersion = csVersion[-2].split(' ')[-2]

    section_dict = {sections[i]: sections[i + 1] for i in range(0, len(sections), 2)}

    # Dictionaries to hold dataframes
    svdf = pd.DataFrame()
    vardf = pd.DataFrame()
    infodf = pd.DataFrame()
    qcdf = pd.DataFrame()

    qcinfo = {}
    qcinfo['MAPPING/ALIGNING SUMMARY: Total bases'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Total input reads'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Insert length: mean'] = 'NA'
    qcinfo['COVERAGE SUMMARY: Average alignment coverage over genome'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Number of duplicate marked reads (%)'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Mapped reads (%)'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Mismatched bases R1 (%)'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Mismatched bases R2 (%)'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Q30 bases R1 (%)'] = 'NA'
    qcinfo['MAPPING/ALIGNING SUMMARY: Q30 bases R1 (%)'] = 'NA'

    ploidy = {'Purity':'NA','Ploidy':'NA'}
    caseInfo = {'caseid':report_id,'version':csVersion,'date':date}
    for k in ['mrn', 'accession', 'specimen_type', 'dob', 'sex', 'runid','instrument', 'flowcell']:
        caseInfo[k] = 'NA'

    for section_name, section_content in section_dict.items():
        # Split the section by lines
        lines = section_content.strip().split("\n")
        
        # Extract case information
        if section_name == 'CHROMOSEQ CASE INFORMATION':
            for line in lines:
                key, value = line.split("\t")[0:2]
                key = key.strip(':').lower().replace(' ','_')
                caseInfo[key] = value            

        # Extract data for different sections
        elif section_name in ['COPY NUMBER ALTERATIONS', 'RECURRENT TRANSLOCATIONS', 'OTHER STRUCTURAL VARIANTS', 'NOVEL/FILTERED SV AND CNV VARIANTS']:
            df = pd.DataFrame()

            data = [line.strip().split("\t") for line in lines[1:]]
            # if its the version 2.0 style:
            if 'csyntax' in lines[0]:
                columns = 'type chrom1 pos1 chrom2 pos2 length csyntax psyntax genes filters id abundance info dblookup'.split(' ')
        
                # check to make sure the data and headers line
                if len(data) > 0 and len(data[0]) != len(columns):
                    return(None)                
                df = pd.DataFrame(data, columns=columns)

            else:
                # old style of report
                columns = ['type', 'chrom1', 'pos1', 'chrom2', 'pos2', 'length', 'BANDS', 'knowngenes', 'psyntax', 'csyntax', 'totalgenes', 'filters', 'id', 'abundance', 'info']

                if len(data) > 0 and len(data[0]) != len(columns):
                    return(None)

                df = pd.DataFrame(data, columns=columns)
                df['genes'] = df['knowngenes']
                df.loc[df['type']=='BND','genes'] = df.loc[df['type']=='BND','totalgenes']

            df['caseid'] = report_id
            df['date'] = date
            df['version'] = csVersion
            df['category'] = categories[section_name]
            df = df[['caseid', 'date', 'version', 'category','type','chrom1','pos1','chrom2','pos2','length','csyntax','psyntax','genes','filters','id','abundance','info']]
            svdf = pd.concat([svdf,df],axis=0).reindex()
        
        elif section_name in ['GENE MUTATIONS', 'FILTERED GENE MUTATIONS']:
            # old style of report
            columns = ['type', 'chrom', 'pos', 'ref', 'alt', 'gene', 'consequence', 'csyntax', 'psyntax', 'exon', 'filters', 'id', 'vaf', 'altreads', 'coverage', 'popaf']
            # v2.0 style
            if 'csyntax' in lines[0]:
                columns = 'type filters chrom pos ref alt gene transcript consequence csyntax psyntax exon popaf annotations coverage altreads vaf dblookup'.split(' ')

            data = [line.split("\t") for line in lines[1:]]

            if len(data) > 0 and len(data[0]) != len(columns):
                return(None)

            df = pd.DataFrame(data, columns=columns)
            df['caseid'] = report_id
            df['date'] = date
            df['version'] = csVersion
            df['category'] = categories[section_name]
            if 'transcript' not in df.columns:
                df['transcript'] = 'NA'
            if 'annotations' not in df.columns:
                df['annotations'] = 'NA'

            df = df[['caseid','date', 'version','category','type','filters', 'chrom', 'pos', 'ref','alt','gene','transcript','consequence','csyntax','psyntax','exon','popaf','annotations', 'coverage','altreads','vaf']]
            vardf = pd.concat([vardf,df],axis=0).reindex()

        elif section_name in ['SEQUENCING QC','CHROMOSEQ QC']:
            for line in lines:
                key, value = line.split("\t")[0:2]
                key = key.strip(':')
                value = value.replace('%','')
                qcinfo[key] = value
        
            # if old style of report
            if 'TOTAL READS' in qcinfo.keys():
                if 'B' in qcinfo['TOTAL READS'] or float(qcinfo['TOTAL READS']) < 1e6:
                    qcinfo['MAPPING/ALIGNING SUMMARY: Total input reads'] = float(qcinfo['TOTAL READS'][:-1]) * 1e9
                else:
                    qcinfo['MAPPING/ALIGNING SUMMARY: Total input reads'] = float(qcinfo['TOTAL READS'][:-1])                    

            if 'TOTAL GIGABASES' in qcinfo.keys():
                qcinfo['MAPPING/ALIGNING SUMMARY: Total bases'] = float(qcinfo['TOTAL GIGABASES']) * 1e9

            if 'MEAN INSERT SIZE' in qcinfo.keys():
                qcinfo['MAPPING/ALIGNING SUMMARY: Insert length: mean'] = float(qcinfo['MEAN INSERT SIZE'])
            
            if 'AVERAGE COVERAGE' in qcinfo.keys():
                qcinfo['COVERAGE SUMMARY: Average alignment coverage over genome'] = float(qcinfo['MEAN INSERT SIZE'])

            if 'DUPLICATES' in qcinfo.keys():
                qcinfo['MAPPING/ALIGNING SUMMARY: Number of duplicate marked reads (%)'] = qcinfo['DUPLICATES']

            if 'MAPPED READS' in qcinfo.keys():
                qcinfo['MAPPING/ALIGNING SUMMARY: Number of duplicate marked reads (%)'] = qcinfo['MAPPED READS']

            df = pd.DataFrame.from_dict([qcinfo])
            df['caseid'] = report_id
            df['date'] = date
            df['version'] = csVersion
            df = df[['caseid','date','version','MAPPING/ALIGNING SUMMARY: Total input reads',
                         'MAPPING/ALIGNING SUMMARY: Total bases',
                         'MAPPING/ALIGNING SUMMARY: Insert length: mean',
                         'COVERAGE SUMMARY: Average alignment coverage over genome',
                         'MAPPING/ALIGNING SUMMARY: Number of duplicate marked reads (%)',
                         'MAPPING/ALIGNING SUMMARY: Mapped reads (%)',
                         'MAPPING/ALIGNING SUMMARY: Mismatched bases R1 (%)',
                         'MAPPING/ALIGNING SUMMARY: Mismatched bases R2 (%)',
                         'MAPPING/ALIGNING SUMMARY: Q30 bases R1 (%)',
                         'MAPPING/ALIGNING SUMMARY: Q30 bases R1 (%)']]
            df.columns = 'caseid date version total_reads total_bases insert_size coverage duplicate_rate mapped_reads mismatched_bases_R1 mismatched_bases_R2 q30_bases_R1 q30_bases_R2'.split(' ')
            qcdf = pd.concat([qcdf,df],axis=1)

        elif section_name == 'Haplotect Contamination Estimate':
            names = [ 'Haplotect.'+ x for x in lines[0].strip().split('\t')[1:] ]
            values = lines[1].split('\t')[1:]
            df = pd.DataFrame.from_dict([dict(zip(names,values))])
            df.columns = 'haplotect_sites haplotect_informative_sites haplotect_coverage haplotect_total_reads haplotect_contaminating_reads haplotect_contamination_fraction haplotect_contamination_estimate haplotect_contamination_confidence'.split(' ')
            qcdf = pd.concat([qcdf,df],axis=1)
            
        elif section_name == 'PLOIDY AND PURITY':
            for line in lines:
                key, value = line.split("\t")[0:2]
                key = key.strip(':')
                ploidy[key] = value

    if 'dob' not in caseInfo.keys():
        caseInfo['dob'] = '1/1/2000'
    infodf = pd.DataFrame.from_dict([caseInfo])
    if all(column in infodf.columns for column in ['mrn', 'accession', 'specimen_type', 'dob', 'sex', 'runid','instrument', 'flowcell']) is False:
        return(None)
            
    infodf = infodf[['caseid','date','version','mrn', 'accession', 'specimen_type', 'dob', 'sex', 'runid','instrument', 'flowcell']]

    qcdf['ploidy'] = ploidy['Ploidy']
    qcdf['purity'] = ploidy['Purity']
    
    svdf['abundance'] = svdf['abundance'].str.replace(r'%|>','',regex=True).astype('float')
    vardf['vaf'] = vardf['vaf'].str.rstrip('%').astype('float')
    vardf['popaf'] = vardf['popaf'].str.rstrip('%').astype('float',errors='ignore')

    dfs = {'variants':vardf,'svs':svdf,'qc':qcdf,'cases':infodf}
    
    return dfs


def parse_cs_report_json(filename):
    if filename:
        print('it works')



def main():

    parser = argparse.ArgumentParser(description='ChromoSeq assay database utilities.')
    subparsers = parser.add_subparsers(title='commands', dest='command')

    # Add command
    add_parser = subparsers.add_parser('add', help='Add ChromoSeq results to database')
    add_parser.add_argument('database', help='ChromoSeq database')
    add_parser.add_argument('reports', nargs='*', help='ChromoSeq reports')
    add_parser.add_argument('-f', '--format', choices=['text', 'json'], default='text', help='Input file(s) are json')
    add_parser.add_argument('-b','--backup', default=False,action='store_true',help='Backup database if >30 days since last backup.')
    add_parser.add_argument('-F','--forcebackup', default=False,action='store_true',help='Force backup regardless of last backup.')

    # Parse command
    parse_parser = subparsers.add_parser('parse', help='Parse ChromoSeq results and print to a file')
    parse_parser.add_argument('reports', nargs='*', help='ChromoSeq reports')
    parse_parser.add_argument('-t', '--table', required=True, choices=['cases', 'qc', 'svs', 'variants'], default='cases', help='Section of the report to parse')
    parse_parser.add_argument('-f', '--format', choices=['text', 'json'], default='text', help='Input file(s) are json')
    parse_parser.add_argument('-o','--outfile', required=False,help='Outfile for parsed output')

    # Dump command
    dump_parser = subparsers.add_parser('dump', help='Dump ChromoSeq database table.')
    dump_parser.add_argument('database', help='ChromoSeq database')
    dump_parser.add_argument('-t','--table', required=True, choices=['cases','qc','svs','variants','history'],help='ChromoSeq table name')
    dump_parser.add_argument('-o','--outfile', required=False,help='Outfile for parsed output')

    # Backup command
    backup_parser = subparsers.add_parser('backup', help='Dump ChromoSeq database table.')
    backup_parser.add_argument('database', help='ChromoSeq database')

    args = parser.parse_args()

    # Check if the database file exists and create it if not
    if args.command == 'add':
        if not os.path.exists(args.database):
            create_db(args.database)

        # Connect to the database        
        engine = create_engine(f'sqlite:///{args.database}')

        # Backup if needed
        if args.backup is True:
            backup_if_needed(args.database, engine, args.forcebackup)

        for file in args.reports:

            dfs = {}
            if args.format == 'text':
                dfs = parse_cs_report_text(file)

            elif args.format == 'json':
                dfs = parse_cs_report_json(file)

            if dfs is None:
                print(f'{file} could not be parsed',file=sys.stderr)
                continue

            if dfs['cases'].shape[0] > 0:
                add_records_from_dataframe(dfs['cases'], cases, engine)

            if dfs['svs'].shape[0] > 0:
                add_records_from_dataframe(dfs['svs'], svs, engine)

            if dfs['variants'].shape[0] > 0:
                add_records_from_dataframe(dfs['variants'], variants, engine)

            if dfs['qc'].shape[0] > 0:
                add_records_from_dataframe(dfs['qc'], qc, engine)

    elif args.command == 'dump':
        if not os.path.exists(args.database):
            raise ValueError(f'Database file {args.database} does not exist')
                # Connect to the database

        engine = create_engine(f'sqlite:///{args.database}')

        df = pd.DataFrame()
        if args.table == 'cases':
            df = query_table(cases, engine)

        elif args.table == 'qc':
            df = query_table(qc, engine)

        elif args.table == 'svs':
            df = query_table(svs, engine)

        elif args.table == 'variants':
            df = query_table(variants, engine)

        elif args.table == 'history':
            df = query_table(history, engine)

        df.to_csv(sys.stdout,sep='\t',index=False,header=True)

    elif args.command == 'parse':

        df = pd.DataFrame()        
        for file in args.reports:
            dfs = {}
            if args.format == 'text':
                dfs = parse_cs_report_text(file)

            elif args.format == 'json':
                dfs = parse_cs_report_json(file)

            if dfs is None:
                continue

            df = pd.concat([df,dfs[args.table]],axis=0,ignore_index=True)

        if args.outfile is not None:
            df.to_csv(args.outfile,sep='\t',index=False,header=True)
        else:
            df.to_csv(sys.stdout,sep='\t',index=False,header=True)
        
    elif args.command == 'backup':
        engine = create_engine(f'sqlite:///{args.database}')
        backup_if_needed(args.database, engine, True)


if __name__ == "__main__":
    main()

