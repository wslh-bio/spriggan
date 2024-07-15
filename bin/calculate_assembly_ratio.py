#!/usr/bin/env python3
import argparse
import os
import logging
import sys


logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

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
    parser.add_argument('-s', '--sample_name',
        metavar='sample_name', 
        type=str, 
        help='Sample name', 
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

def does_quast_exist(quast_report):
    # Check if quast report exists
    if not os.path.isfile(quast_report):
        logging.CRITICAL(f"No quast report '{quast_report}' found, cannot continue")

        sys.exit(1)

def initialize_variables():

    # Initialize variables
    taxid = "NA"
    stdev = "NA"
    stdevs = "NA"
    assembly_length = "NA"
    expected_length = "NA"
    total_tax = "NA"

    return taxid, stdev, stdevs, assembly_length, expected_length, total_tax

def process_database_paths(path_database, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    # Process database path
    NCBI_ratio_date = path_database[::-1].split('_')[0][::-1].split('.')[1]

    if os.path.isfile(path_database):

        db_path_update = path_database + "_update.txt"

        with open(path_database, 'r') as infile, open(db_path_update, 'w') as outfile:
            for line in infile:
                outfile.write(line.capitalize().replace('[', '').replace(']', ''))

        NCBI_ratio = db_path_update

        return NCBI_ratio, NCBI_ratio_date

    else:

        logging.CRITICAL("No ratio DB, exiting")
        logging.DEBUG("Writing NA for all information in output files.")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -2")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write("Tax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(1)

def check_quast_stats(quast_report, NCBI, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    # Check if quast assembly stats exist
    with open(quast_report, 'r') as infile:

        lines = infile.readlines()

        if len(lines) >= 17:

            assembly_length = lines[16].split('\t')[2]
            sample_gc_percent = lines[17].split('\t')[2]

            return assembly_length, sample_gc_percent

        else:

            logging.CRITICAL("No quast exists, cannot continue")

            with open(f"{sample_name}_Assembly_ratio_{NCBI}.txt", 'w') as outfile:
                outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -2")

            with open(f"{sample_name}_GC_content_{NCBI}.txt", 'w') as outfile:
                outfile.write("Tax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

            sys.exit(1)

def does_tax_file_exist(tax):

    # Check for tax summary file
    if not os.path.isfile(tax):

        logging.CRITICAL("No Tax file to find accession for lookup, exiting")

        sys.exit(1)

def process_NCBI_and_tax(taxonomy_to_compare, tax):    # Initialize variables for tax lookup

    genus = None
    species = None

    # Process taxonomy information
    if not taxonomy_to_compare:

        with open(tax, 'r') as infile:
            lines = infile.readlines()

            if len(lines) >= 8:

                genus = lines[6].split('\t')[1].strip()

                if genus == "":
                    genus = "No genus found"

                species = lines[7].split('\t')[1].strip()

                if "sp." in species:
                    species = species.replace("sp. ", "sp.").capitalize().replace(" ", "-")

                if species == "":
                    species = "No species found"

                total_tax = f"{genus} {species}"

                return total_tax, genus, species

    else:

        in_genus = taxonomy_to_compare.split()[0].capitalize()
        in_species = taxonomy_to_compare.split()[1].lower().capitalize()
        genus = in_genus
        species = in_species
        total_tax = f"{genus} {species}    (selected manually)"

    # Initialize variables for matching
    found = False

    return total_tax, genus, species, found

def search_ncbi_ratio_file(NCBI_ratio, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax):

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

                if reference_count < 10:
                    stdev = "Not calculated on species with n<10 references"
                    stdevs = "NA"

                else:

                    if int(assembly_length) > expected_length:
                        bigger = int(assembly_length)
                        smaller = expected_length

                    else:
                        smaller = int(assembly_length)
                        bigger = expected_length

                    stdevs = (bigger - smaller) / stdev

                gc_min = line[7]
                gc_max = line[8]
                gc_mean = line[10]
                gc_count = line[12]

                if gc_count < 10:
                    gc_stdev = "Not calculated on species with n<10 references"

                else:
                    gc_stdev = line[11]

                with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
                    outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

                found = True

                return stdev, gc_stdev, gc_min, gc_max, gc_mean, gc_count

    # Handle unmatched cases
    if not found:

        logging.INFO(f"No match found for '{genus} {species}'")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -1")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(0)

def calculate_ratio(sample_name, NCBI_ratio_date, expected_length, total_tax, taxid, assembly_length):

    # Calculate and print ratio
    if expected_length == "NA" or not expected_length:

        logging.INFO("No expected length was found to compare to")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: NA\nRatio: -1")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(0)

    elif assembly_length == "NA" or not assembly_length:

        logging.INFO("No assembly length was found to compare with")

        with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: NA\nActual_length: NA\nExpected_length: {expected_length}\nRatio: -2")

        with open(f"{sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: NA")
        sys.exit(0)

    ratio = float(assembly_length) / float(expected_length)

    logging.INFO(f"Actual - {assembly_length}")
    logging.INFO(f"Expected - {expected_length}")
    logging.INFO(f"Ratio - {ratio:.4f}")

    return ratio

def write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, stdevs, assembly_length, expected_length, ratio):

    with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_St.Dev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: {ratio:.4f}")

def print_version(version):

    if version:
        logging.INFO("calculate_assembly_ratio.py: 2.0")
        sys.exit(0)

def main(args=None):
    args = parse_args(args)

    if args.version == True:
        print_version(args.version)

    does_quast_exist(args.quast_report)
    taxid, stdev, stdevs, assembly_length, expected_length, total_tax = initialize_variables()
    NCBI_ratio, NCBI_ratio_date = process_database_paths(args.path_database, args.sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax)
    does_tax_file_exist(args.tax_file)
    total_tax, genus, species = process_NCBI_and_tax(args.taxonomy_to_compare, args.tax_file)
    ratio = calculate_ratio(args.sample_name, NCBI_ratio, expected_length, total_tax, taxid, assembly_length)
    write_output(args.sample_name, NCBI_ratio_date, total_tax, taxid, stdev, stdevs, assembly_length, expected_length, ratio)

if __name__ == "__main__":
    sys.exit(main())