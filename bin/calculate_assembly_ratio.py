#!/usr/bin/env python3
import argparse
import re
import os
import logging
import sys
import pandas as pd
import urllib.request
import numpy as np
import statistics
from datetime import datetime


logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

DEFAULT_REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

timestamp = datetime.now().strftime("%Y%m%d")

def parse_args(args=None):
    Description='Compare local assembly to expected assembly size based on taxonomy.'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-d', '--path_database',
        metavar='path_to_database_file', 
        type=str,
        default=DEFAULT_REFSEQ_URL,
        help='FTP path of the Refseq assembly summary file (default: {DEFAULT_REFSEQ_URL})', 
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
        help='Tax file from Kraken', 
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
    taxid = "NA" #DATABASE
    stdev = "NA" #CALCULATED
    stdevs = "NA" #CALCULATED
    assembly_length = "NA" #spriggan report, quast_results, 
    expected_length = "NA" #CALCULATED
    total_tax = "NA" #DATABASE

    return taxid, stdev, stdevs, assembly_length, expected_length, total_tax


def is_float(value):
    """
    Helper function to handle strings in genome_size column of refseq assembly file.
    Returns True if the value can be converted to a float.
    """
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False
    

def parse_gc_percent(value):
    """
    Some fields in the gc_percent column of the refseq assembly file contain values over 100%.
    Helper function that gets rid of values over 100%.
    """
    try:
        val = float(value)
        if val < 0:  # extremely unlikely scenario but handling anyways
            return None
        elif val > 100:  # Not sure if this is the best way to handle this but it's the cleanest
            return None
        else:
            return val
    except (ValueError, TypeError):
        return None


def compute_taxid_genome_stats(url, target_taxid, sample_name, assembly_length, total_tax, sample_gc_percent):
    """
    Stream through the NCBI assembly_summary_refseq.txt file,
    compute mean genome_size (after IQR filtering) and mean gc_percent
    for all records matching the target_taxid.
    """

    logging.info(f"Fetching NCBI assembly summary for taxid {target_taxid} ...")

    genome_sizes = []
    gc_percents = []

    try:
        with urllib.request.urlopen(url) as response:
            for raw_line in response:
                line = raw_line.decode("utf-8").strip()
                if not line or line.startswith("#"):
                    continue

                cols = line.split("\t")
                if len(cols) <= 27:  # At minimum, column count should go to column 27 (zero-based numbering), which is 'gc_percent'. 
                    continue

                taxid = cols[5].strip()
                if taxid != str(target_taxid):
                    continue

                genome_size_val = cols[25].strip()
                gc_percent_val = parse_gc_percent(cols[27].strip())  # Remove all gc_percent values over 100

                # Handle strings in genome_size and validate numeric fields
                if is_float(genome_size_val) and gc_percent_val is not None:
                    genome_size = float(genome_size_val)
                    gc_percent = float(gc_percent_val)
                    genome_sizes.append(genome_size)
                    gc_percents.append(gc_percent)
                else:
                    logging.debug(
                        f"Skipping malformed entry for taxid {taxid}: genome_size='{genome_size_val}', gc_percent='{gc_percent_val}'"
                    )

        if not genome_sizes:
            logging.warning(f"No valid genome entries found for taxid {target_taxid}")
            return None, None

        # --- IQR filtering on genome_size and gc_percent ---
        Q1 = np.percentile(genome_sizes, 25)
        Q3 = np.percentile(genome_sizes, 75)
        IQR = Q3 - Q1

        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        filtered_sizes = [x for x in genome_sizes if lower_bound <= x <= upper_bound]  # Remove outliers

        Q1_gc = np.percentile(gc_percents, 25)
        Q3_gc = np.percentile(gc_percents, 75)
        IQR_gc = Q3_gc - Q1_gc

        gc_lower_bound = Q1_gc - 1.5 * IQR_gc
        gc_upper_bound = Q3_gc + 1.5 * IQR_gc

        filtered_gc = [x for x in gc_percents if gc_lower_bound <= x <= gc_upper_bound]  # Remove GC% outliers

        if not filtered_sizes or not filtered_gc:
            logging.warning(f"All values filtered out for taxid {target_taxid}")
            return None

        # Final stats
        expected_length = statistics.mean(filtered_sizes)  # Mean genome size
        # For continuity, only calculate std dev for species with 10 or more references
        if len(filtered_sizes) >= 10:
            stdev_genome_size = statistics.stdev(filtered_sizes)
            if int(assembly_length) > int(expected_length):
                bigger = int(assembly_length)
                smaller = int(expected_length)
            else:
                smaller = int(assembly_length)
                bigger = int(expected_length)

            logging.debug("Calculating the standard deviations")
            stdevs = (bigger - smaller) / stdev_genome_size

        else:
            stdev_genome_size = "Not calculated on species with n<10 references"
            stdevs = "NA"
        
        ## GC
        species_gc_mean = statistics.mean(gc_percents)
        gc_min = min(gc_percents)
        gc_max = max(gc_percents)
        gc_count = len(gc_percents)
        species_gc_percent_stdev = statistics.stdev(filtered_gc) if len(filtered_gc) >= 10 else "Not calculated on species with n<10 references" 


        logging.info(
            f"Taxid {target_taxid}: mean genome_size (IQR-filtered) = {expected_length:.2f}, "
            f"mean GC% = {species_gc_mean:.2f} (n={len(filtered_sizes)})"
        )

        with open(f"{sample_name}_GC_content_{timestamp}.txt", 'w') as outfile:
                    outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {species_gc_percent_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {species_gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

        return stdev_genome_size, species_gc_percent_stdev, gc_min, gc_max, species_gc_mean, gc_count, stdevs, expected_length, taxid
    
    except Exception as e:
        logging.error(f"Error computing taxid genome stats for {target_taxid}: {e}")
        return None


def extract_sample_name(quast_report):
    logging.debug("Extracting sample name from quast report")

    sample_name = quast_report.split('.')[0]
    sample_name = sample_name.split('/')[-1]

    return sample_name

def handle_missing_database_paths(path_database, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    logging.debug("Function to handle erroneous database path passed into main function")
    if not os.path.isfile(path_database):
        logging.critical("No ratio database found, exiting")
        logging.debug("Writing NA for all information in output files.")

        with open(f"{sample_name}_Assembly_ratio_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: -2\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(1)

def check_quast_stats(quast_report, NCBI, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax):

    logging.debug("Checking quast results.")

    if quast_report:

        # Check if quast assembly stats exist
        with open(quast_report, 'r') as infile:

            logging.debug("Setting up column names here")
            assembly_length_row_name = "Total length (>= 0 bp)"
            gc_percent_row_name = "GC (%)"

            logging.debug("Going through quast file to get result for assembly length and sample gc")
            df = pd.read_csv(infile,sep = '\t')
            assembly_length = df[assembly_length_row_name][0]
            sample_gc_percent = df[gc_percent_row_name][0]

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

    logging.debug("Initializing blank variables")
    genus = None
    species = None

    logging.debug("If user did not enter a taxonomy to compare to")
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
        logging.debug("If user has a taxonomy they want to compare this to")
        in_genus = taxonomy_to_compare.split()[0].capitalize()
        in_species = taxonomy_to_compare.split()[1].lower().capitalize()
        genus = in_genus
        species = in_species

        total_tax = f"{genus} {species}    (selected manually)"

        # Initialize variables for matching
        found = False

        return total_tax, genus, species, found


def calculate_ratio(sample_name, expected_length, total_tax, taxid, assembly_length, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdev):

    if expected_length == "NA" or not expected_length:

        logging.info("No expected length was found to compare to")

        with open(f"{sample_name}_Assembly_ratio_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: NA\nRatio Actual:Expected: -1\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")

        sys.exit(0)

    elif assembly_length == "NA" or not assembly_length:

        logging.info("No assembly length was found to compare with")

        with open(f"{sample_name}_Assembly_ratio_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: NA\nActual_length: NA\nExpected_length: {expected_length}\nRatio Actual:Expected: -2\nRatio Expected:Actual: NA")

        with open(f"{sample_name}_GC_content_{timestamp}.txt", 'w') as outfile:
            outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: NA")
        sys.exit(0)

    logging.debug("Calculating the assembly and expected ratios")
    ratio_a_e = float(assembly_length) / float(expected_length)
    ratio_e_a = float(expected_length) / float(assembly_length)

    logging.info(f"Actual - {assembly_length}")
    logging.info(f"Expected - {expected_length}")
    logging.info(f"Ratio Actual:Expected - {ratio_a_e:.4f}")
    logging.info(f"Ratio Expected:Actual - {ratio_e_a:.4f}")

    return ratio_a_e, ratio_e_a

def write_output(sample_name, total_tax, taxid, stdev, stdevs, assembly_length, expected_length, ratio_a_e, ratio_e_a):

    logging.debug("Writing the output file")
    with open(f"{sample_name}_Assembly_ratio_{timestamp}.txt", 'w') as outfile:
        outfile.write(f"Sample: {sample_name}\nTax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_St.Dev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio Actual:Expected: {ratio_a_e}\nRatio Expected:Actual: {ratio_e_a}")

def print_version(version):
    logging.debug("Took this version from the original script")
    if version:
        logging.info("calculate_assembly_ratio.py: 2.0")

def main(args=None):
    args = parse_args(args)

    if args.version:
        print_version(args.version)

    #Initializing the variables
    taxid, stdev, stdevs, assembly_length, expected_length, total_tax = initialize_variables()

    #Extracting sample name
    sample_name = extract_sample_name(args.quast_report)

    #Grabbing assembly length and gc percentage from quast file
    assembly_length, sample_gc_percent = check_quast_stats(args.quast_report, args.tax_file, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax)

    #Getting database names and dates
    handle_missing_database_paths(args.path_database, sample_name, taxid, stdev, stdevs, assembly_length, expected_length, total_tax)

    #Getting taxonomy info
    total_tax, genus, species, found = process_NCBI_and_tax(args.taxonomy_to_compare, args.tax_file, sample_name)

    #Grabbing stats 
    # stdev, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdevs, expected_length, taxid = search_ncbi_ratio_file(NCBI_ratio_file, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax, sample_gc_percent, found)

    #TODO Replace with compute_taxid_genome_stats()
    # result = search_ncbi_ratio_file(
    #     NCBI_ratio_file, genus, species, assembly_length,
    #     sample_name, NCBI_ratio_date, total_tax,
    #     sample_gc_percent, found
    # )

    if result is None:
        logging.warning(
            f"No NCBI assembly stats found for '{genus} {species}'in {args.path_database}; proceeding with default values for sample '{sample_name}'."
        )
        # Assign placeholder values so output can still be written
        taxid = "None"
        stdev = "None"
        stdevs = "None"
        assembly_length_str = str(assembly_length) if assembly_length else "None"
        expected_length = "None"
        total_tax = "None"
        ratio_a_e = "None"
        ratio_e_a = "None"

        # Write placeholder Assembly ratio file
        with open(f"{sample_name}_Assembly_ratio_{timestamp}.txt", 'w') as outfile:
            outfile.write(
                f"Sample: {sample_name}\n"
                f"Tax: {total_tax}\n"
                f"NCBI_TAXID: {taxid}\n"
                f"Species_St.Dev: {stdev}\n"
                f"Isolate_St.Devs: {stdevs}\n"
                f"Actual_length: {assembly_length_str}\n"
                f"Expected_length: {expected_length}\n"
                f"Ratio Actual:Expected: {ratio_a_e}\n"
                f"Ratio Expected:Actual: {ratio_e_a}\n"
            )

        # Write placeholder GC content file
        with open(f"{sample_name}_GC_content_{timestamp}.txt", 'w') as outfile:
            outfile.write(
                f"Sample: {sample_name}\n"
                f"Tax: {total_tax}\n"
                f"NCBI_TAXID: {taxid}\n"
                f"Species_GC_StDev: None\n"
                f"Species_GC_Min: None\n"
                f"Species_GC_Max: None\n"
                f"Species_GC_Mean: None\n"
                f"Species_GC_Count: None\n"
                f"Sample_GC_Percent: None\n"
            )

    else:
        stdev, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdevs, expected_length, taxid = result

        # Calculate ratios
        ratio_a_e, ratio_e_a = calculate_ratio(
            sample_name, timestamp, expected_length, total_tax, taxid,
            assembly_length, gc_stdev, gc_min, gc_max, gc_mean, gc_count, stdev
        )

        # Write final output files
        write_output(
            sample_name, timestamp, total_tax, taxid, stdev, stdevs,
            assembly_length, expected_length, ratio_a_e, ratio_e_a
            )
    logging.info("Finished writing assembly ratio file and GC content file.")

if __name__ == "__main__":
    sys.exit(main())
