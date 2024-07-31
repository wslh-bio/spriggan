#!/usr/bin/env python3
import argparse
import re
import os
import logging
import sys
import pandas as pd


logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='Compare local assembly to expected assembly size based on taxonomy.'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-d', '--path_database',
        metavar='path_to_database_file', 
        type=str, 
        help='Path to sorted database with NCBI entries statistics', 
        required=True
        )
    parser.add_argument('-q', '--quast_report',
        metavar='report.tsv', 
        type=str, 
        help='Quast report.tsv file', 
        required=True
        )
    parser.add_argument('-t', '--tax_file',
        metavar='tax_file', 
        type=str, 
        help='Tax file from determine_taxaID.sh output', 
        required=True
        )
    parser.add_argument('-f', '--taxonomy_to_compare',
        metavar='"genus species"', 
        type=str, 
        help='Specific taxonomy to compare against in the database', 
        default=None
        )
    parser.add_argument('-V', '--version',
        action='store_true', 
        help='Print version and exit'
        )
    return parser.parse_args(args)

def initialize_variables():
    logging.debug("Initializing variables")

    # Initialize variables
    taxid = "NA"
    stdev = "NA" #CALCULATED
    stdevs = "NA" #CALCULATED
    assembly_length = "NA" #spriggan report, quast_results, 
    expected_length = "NA" #CALCULATED
    total_tax = "NA"

    return taxid, stdev, stdevs, assembly_length, expected_length, total_tax

def extract_sample_name(quast_report):
    logging.debug("Extracting sample name from quast report")
    sample_name = quast_report.split('.')[0]
    sample_name = sample_name.split('/')[-1]

    return sample_name

def process_database_paths(path_database, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    logging.debug("Processing database dates and paths.")

    # Process database path
    get_date = path_database
    match = re.search(r"(\d{8})", get_date)
    NCBI_ratio_date = match.group(1)

    file_name_txt = os.path.basename(path_database)
    file_name = file_name_txt.strip(".txt")

    dir_name = os.path.dirname(path_database) + "/"

    if os.path.isfile(path_database):

        db_path_update = dir_name + file_name + "_update.txt"

        with open(path_database, 'r') as infile, open(db_path_update, 'w') as outfile:
            for line in infile:
                outfile.write(line.capitalize().replace('[', '').replace(']', ''))

        NCBI_ratio = db_path_update

        return NCBI_ratio, NCBI_ratio_date

    else:

        logging.critical("No ratio DB, exiting")
        logging.debug("Writing NA for all information in output files.")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: -2\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(1)

def check_quast_stats(quast_report, NCBI, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    logging.debug("Checking quast results.")

    if quast_report:

        # Check if quast assembly stats exist
        with open(quast_report, 'r') as infile:

            df = pd.read_csv(infile,sep='\t')

            assembly_length = df.iloc[6].values[1]
            sample_gc_percent = df.iloc[15].values[1]


            return assembly_length, sample_gc_percent

    else:

        logging.critical("No quast exists, cannot continue")

        with open(f"{sample_name}_Assembly_ratio_{NCBI}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: -2\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{NCBI}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(1)

def process_NCBI_and_tax(taxonomy_to_compare, tax, sample_name):

    logging.debug("Checking for taxonomy information in taxonomy file.")

    genus = None
    species = None

    # Process taxonomy information
    if not taxonomy_to_compare:

        df = pd.read_csv(tax,sep='\t')

        if df['Sample'].str.contains(sample_name).any():

            found = True

            result = df.loc[df['Sample'] == sample_name,'Primary Species (%)'].values[0]

            genus = result.split(' ')[0]

            if genus == "":
                genus = "No genus found"

            species = result.split(' ')[1]

            if 'sp.' in species:
                species = species.replace('sp. ', 'sp.')

            if species == "":
                species = "No species found"

            total_tax = f"{genus} {species}"

            return total_tax, genus, species, found

    else:

        in_genus = taxonomy_to_compare.split()[0].capitalize()
        in_species = taxonomy_to_compare.split()[1].lower().capitalize()
        genus = in_genus
        species = in_species

        total_tax = f"{genus} {species}    (selected manually)"

        # Initialize variables for matching
        found = False

        return total_tax, genus, species, found

def search_ncbi_ratio_file(NCBI_ratio, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax, sample_gc_percent, found):

    # Search in NCBI_ratio file
    with open(NCBI_ratio, 'r') as infile:
        for line in infile:

            line = line.strip().split('\t')

            if f"{genus.lower()} {species.lower()}" == line[0].lower():
                taxid = line[19]

                if taxid == -2:
                    taxid = "No mode available when determining tax id"

                elif taxid == -1:
                    taxid = "No tax id given or empty when making lookup"

                expected_length = int(1000000 * float(line[4])) // 1
                reference_count = line[6]
                stdev = int(1000000 * float(line[5])) // 1

                if int(reference_count) < 10:
                    stdev = "Not calculated on species with n<10 references"
                    stdevs = "NA"

                else:
                    if int(assembly_length) > int(expected_length):
                        bigger = int(assembly_length)
                        smaller = int(expected_length)

                    else:
                        smaller = int(assembly_length)
                        bigger = int(expected_length)

                    stdevs = (bigger - smaller) / stdev

                gc_min = line[7]
                gc_max = line[8]
                gc_mean = line[10]
                gc_count = line[12]

                if int(gc_count) < 10:
                    gc_stdev = "Not calculated on species with n<10 references"

                else:
                    gc_stdev = line[11]

                with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
                    outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

                found = True

                return stdev, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdevs, expected_length, taxid

    # Handle unmatched cases
    if not found:

        logging.info(f"No match found for '{genus} {species}'")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: -1\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(0)

def calculate_ratio(sample_name, NCBI_ratio_date, expected_length, total_tax, taxid, assembly_length, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdev):

    # Calculate and print ratio
    if expected_length == "NA" or not expected_length:

        logging.info("No expected length was found to compare to")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: NA\nRatio Actual:Expected: -1\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(0)

    elif assembly_length == "NA" or not assembly_length:

        logging.info("No assembly length was found to compare with")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: NA\nActual_length: NA\nExpected_length: {expected_length}\nRatio Actual:Expected: -2\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: NA")
        sys.exit(0)

    ratio_a_e = float(assembly_length) / float(expected_length)
    ratio_e_a = float(expected_length) / float(assembly_length)

    logging.info(f"Actual - {assembly_length}")
    logging.info(f"Expected - {expected_length}")
    logging.info(f"Ratio Actual:Expected - {ratio_a_e:.4f}")
    logging.info(f"Ratio Expected:Actual - {ratio_e_a:.4f}")

    return ratio_a_e, ratio_e_a

def write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, stdevs, assembly_length, expected_length, ratio_a_e, ratio_e_a):

    with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_St.Dev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: {ratio_a_e}\nRatio Expected:Actual: {ratio_e_a}")

def print_version(version):

    if version:
        logging.info("calculate_assembly_ratio.py: 2.0")

def main(args=None):
    args = parse_args(args)

    if args.version == True:
        print_version(args.version)

    #Initializing the variables
    taxid, stdev, stdevs, assembly_length, expected_length, total_tax = initialize_variables()

    #Extracting sample name
    sample_name = extract_sample_name(args.quast_report)

    #Grabbing assembly length and gc percentage from quast file
    assembly_length, sample_gc_percent = check_quast_stats(args.quast_report, args.tax_file, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax)

    #Getting database names and dates
    NCBI_ratio, NCBI_ratio_date = process_database_paths(args.path_database, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax)

    #Getting taxonomy info
    total_tax, genus, species, found = process_NCBI_and_tax(args.taxonomy_to_compare, args.tax_file, sample_name)

    #Grabbing stats 
    stdev, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdevs, expected_length, taxid = search_ncbi_ratio_file(NCBI_ratio, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax, sample_gc_percent, found)

    #Calculating ratio 
    ratio_a_e, ratio_e_a = calculate_ratio(sample_name, NCBI_ratio, expected_length, total_tax, taxid, assembly_length,gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdev)

    #Writing final output
    write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, stdevs, assembly_length, expected_length, ratio_a_e, ratio_e_a)

    logging.info("Finished writing assembly ratio file.")

if __name__ == "__main__":
    sys.exit(main())