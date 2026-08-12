#!/usr/bin/env python3
import sys
import glob
import pandas as pd
from functools import reduce
import argparse
import re
import logging
import csv

def sanitize_sample(sample):
    """
    Strip suffixes if the sample matches internal format (e.g., 25MP######_S01).
    """
    if pd.isna(sample):
        return sample

    s = str(sample)

    # Internal pattern: two digits + MP + six digits
    mp = re.match(r"^\d{2}MP\d{6,}", s)
    if mp:
        return mp.group(0)

    # For non-matching sample names
    return s

def sanitize_primary_species(primary_species):
    """
    Split out the % from primary species and assign to "Primary Species" column
    """
    if pd.isna(primary_species):
        return primary_species

    species = str(primary_species)
    cleaned_species = species.split('(')[0].strip()
    return cleaned_species

def modify_mlst_scheme(row):
    """
    Modify the MLST Scheme column to format MLST<WGS result>_<scheme used>_<organism name>
    """
    scheme_str = row['MLST Scheme']
    species = row['Primary Species (%)'].split(' (')[0]

    # handle missing or empty MLST Scheme
    if pd.isna(scheme_str) or scheme_str.strip() == "":
        return scheme_str

    # split on semicolon if present (for organisms with multiple schemes)
    schemes = scheme_str.split(';') if ';' in scheme_str else [scheme_str]

    # modify the MLST scheme with species appended to each
    modified = [f"{s}_{species}" for s in schemes]

    # join with "_" separator
    return '_'.join(modified)

def create_dataframe(result_files):
    """
    Create the output data frame
    """
    logging.debug(f"Initial result files: {result_files}")

    logging.debug("Getting all tsv files and read them in as data frames")
    tsv_files = glob.glob("*.tsv")

    logging.debug("Setting up df for all result files")
    dfs = []

    logging.debug("Getting files that should not be merged and set them up as ad")
    kraken_ntc_files = glob.glob("*.kraken2.txt")
    db_versions = glob.glob("*.yml")


    for file in tsv_files:
        logging.debug(f"File to be merged: {file}")
        df = pd.read_csv(file, header=0, delimiter='\t')

        logging.debug(f"Sanitizing {file} samples names")
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].apply(sanitize_sample)

        logging.debug(f"Sanitizing {file} primary species")
        if 'Primary Species (%)' in df.columns:
            df['Primary Species'] = df['Primary Species (%)'].apply(sanitize_primary_species)

        dfs.append(df)

    logging.debug("Merging data frames based on sample")
    merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='outer'), dfs)

    logging.debug("Converting sample names to string")
    merged_df['Sample'] = merged_df['Sample'].astype(str)

    logging.debug("Modifying MLST scheme")
    merged_df['MLST Scheme'] = merged_df.apply(modify_mlst_scheme, axis=1)

    return merged_df, kraken_ntc_files, db_versions

def grab_db_versions(db_versions):
    logging.debug("Opening version files to get DB versions")

    for version_file in db_versions:
        with open(version_file, 'r') as programFile:
            for l in programFile.readlines():
                if "kraken DB:" in l.strip():
                    krakenDBVersion = l.strip().split(':')[1].strip()
                if "amrfinder DB:" in l.strip():
                    amrfinderVersion = l.strip().split(':')[1].strip()

    return krakenDBVersion, amrfinderVersion

def kraken_ntc_processing_and_empty_check(kraken_ntc_files, empty_ntcs, merged_df):

    logging.debug("Getting Kraken NTC results")
    if kraken_ntc_files != []:
        kraken_ntc_results = glob.glob("*.kraken2.txt")

        logging.debug("Adding NTC column and calculate Kraken NTC read totals")
        ntc_total_reads = []
        max_ntc_reads = 0

        logging.debug("Reading in kraken NTC files and get # of total reads")
        for file in kraken_ntc_results:
            id = file.split(".kraken2.txt")[0]
            total_reads = 0

            with open(file,'r') as csvfile:
                dialect = csv.Sniffer().sniff(csvfile.read(1024))
                csvfile.seek(0)
                reader = csv.reader(csvfile,dialect)
                for row in reader:
                    if row[3] == "U":
                        total_reads += int(row[1])
                        if int(row[1]) > max_ntc_reads:
                            max_ntc_reads = int(row[1])
                    if "root" in row[5]:
                        total_reads += int(row[1])
                        if int(row[1]) > max_ntc_reads:
                            max_ntc_reads = int(row[1])

        logging.debug("Checking if any NTCs are empty and adds them to the totals.")
        string = ''.join(empty_ntcs)
        filtered = string.strip("[]")
        empty_NTC_list = filtered.split(",")

        for sample in empty_NTC_list:
            if sample != "Empty":
                if sample not in ntc_total_reads:
                    ntc_total_reads.append(f"{sample}: 0")
                    total_reads += 0

        logging.debug("Assigning max reads for ntcs")
        merged_df = merged_df.assign(max_ntc_reads=max_ntc_reads)

        ntc_total_reads.append(f"{id}: {total_reads}")

        logging.debug("Adding NTC totals to data frame")
        merged_df = merged_df.assign(ntc_all_reads=", ".join(ntc_total_reads))

    else:
        logging.debug("Accounting for empty kraken NTC file(s)")
        merged_df = merged_df.assign(ntc_all_reads="999999")
        merged_df = merged_df.assign(max_ntc_reads="999999")

    return merged_df

def assign_versions(merged_df, krakenDBVersion, amrfinderDBVersion, WFVersion):

    logging.debug("Adding kraken DB column")
    merged_df = merged_df.assign(krakenDB=krakenDBVersion)

    logging.debug("Adding AMRFinderPlus DB column")
    merged_df = merged_df.assign(amrfinderDB=amrfinderDBVersion)

    logging.debug("Adding Workflow version column")
    merged_df = merged_df.assign(workflowVersion=WFVersion)

    return merged_df

def assign_run_name(merged_df, WFRunName):

    logging.debug("Using the workflow run name from params for the run column")
    merged_df['Run'] = f"{WFRunName}"

    return merged_df

def rename_columns(merged_df):
    logging.debug("Renaming columns to nicer names")
    merged_df = merged_df.rename(columns={'Contigs':'Contigs (#)',
                                'Average Coverage':'Mean Coverage',
                                'Gene':'AMR','Coverage':'AMR Coverage',
                                'Identity':'AMR Identity',
                                'krakenDB':'Kraken Database Version',
                                'amrfinderDB':'AMRFinderPlus Database Version',
                                'workflowVersion':'Spriggan Version',
                                'ntc_all_reads':'All NTC reads',
                                'max_ntc_reads':'Max NTC read'})
    return merged_df

def reorder_columns(merged_df):
    logging.debug("Putting columns in specific order")
    merged_df = merged_df[['Sample',
                    'Run',
                    'Total Reads',
                    'Reads Removed',
                    'Median Coverage',
                    'Average Coverage',
                    'Contigs',
                    'Assembly Length (bp)',
                    'N50',
                    'Primary Species (%)',
                    'Secondary Species (%)',
                    'Unclassified Reads (%)',
                    'krakenDB',
                    'MLST Scheme',
                    'Gene',
                    'Coverage',
                    'Identity',
                    'Selected AMR Genes',
                    'Selected AMR Genes Coverage',
                    'Selected AMR Genes Identity',
                    'Expected Genome Length',
                    'Genome Length Ratio (Actual/Expected)',
                    'Species GC Content (Mean)',
                    'Sample GC Content (%)',
                    'Primary Species',
                    'ntc_all_reads',
                    'max_ntc_reads',
                    'amrfinderDB',
                    'workflowVersion']]
    return merged_df

def write_output(WFRunName, merged_df):

    logging.info("Writing results to csv file")
    merged_df.to_csv(f'{WFRunName}_spriggan_report.csv', index=False, sep=',', encoding='utf-8')

class CompiledResults(argparse.ArgumentParser):

    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

        sys.exit(1)

if __name__ == "__main__":
    parser = CompiledResults(prog = 'Compiles all SPNTypeID results',
        description='A script to summarize stats',
        epilog='Use with create_report.py --result_files <CH_RESULTS> --workflowRunName <RUN_NAME> --empty_ntc_list <EMPTY_NTC_LIST>'
        )
    parser.add_argument('--result_files',
        nargs="+",
        help='Compiled results from SPNtypeID'
        )
    parser.add_argument('--workflowVersion',
        type=str,
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.'
        )
    parser.add_argument('--workflowRunName',
        type=str,
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.'
        )
    parser.add_argument('--empty_ntc_list',
        nargs="*",
        help='This is determined in the spnetypeid script.'
        )

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    logging.info("Begin compiling all results for final output file.")
    merged_df, kraken_ntc_files, db_versions = create_dataframe(args.result_files)

    krakenDBVersion, amrfinderDBVersion = grab_db_versions(db_versions)

    merged_df = assign_versions(merged_df, krakenDBVersion, amrfinderDBVersion, args.workflowVersion)

    merged_df = kraken_ntc_processing_and_empty_check(kraken_ntc_files, args.empty_ntc_list, merged_df)

    merged_df = assign_run_name(merged_df, args.workflowRunName)

    merged_df = reorder_columns(merged_df)

    merged_df = rename_columns(merged_df)

    write_output(args.workflowRunName, merged_df)
