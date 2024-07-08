#!/usr/bin/env python3
import argparse
import os
import logging


logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

# Parse command line options
parser = argparse.ArgumentParser(
    description='Compare local assembly to expected assembly size based on taxonomy.'
    )
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

args = parser.parse_args()

if args.version:
    logging.INFO("calculate_assembly_ratio.py: 2.0")
    exit(0)

# Check if quast report exists
if not os.path.isfile(args.quast_report):
    logging.CRITICAL(f"No quast report '{args.quast_report}' found, cannot continue")
    exit(1)

# Initialize variables
taxid = "NA"
stdev = "NA"
stdevs = "NA"
assembly_length = "NA"
expected_length = "NA"
total_tax = "NA"

# Process database path
NCBI_ratio_date = args.path_database[::-1].split('_')[0][::-1].split('.')[1]

if os.path.isfile(args.database):

    db_path_update = args.path_database + "_update.txt"

    with open(args.d, 'r') as infile, open(db_path_update, 'w') as outfile:
        for line in infile:
            outfile.write(line.capitalize().replace('[', '').replace(']', ''))

    NCBI_ratio = db_path_update

else:

    logging.CRITICAL("No ratio DB, exiting")

    with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -2")

    with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write("Tax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")
    exit(1)

# Check if quast assembly stats exist
with open(args.quast_report, 'r') as infile:

    lines = infile.readlines()

    if len(lines) >= 17:

        assembly_length = lines[16].split('\t')[2]
        sample_gc_percent = lines[17].split('\t')[2]

    else:

        logging.CRITICAL("No quast exists, cannot continue")

        with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -2")

        with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
            outfile.write("Tax: No genus Found    No species found\nNCBI_TAXID: No Match Found\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")
        exit(1)

# Check for tax summary file
if not os.path.isfile(args.t):
    logging.CRITICAL("No Tax file to find accession for lookup, exiting")
    exit(1)

# Initialize variables for tax lookup
genus = None
species = None

# Process taxonomy information
if not args.taxonomy_to_compare:

    with open(args.t, 'r') as infile:
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

else:

    in_genus = args.taxonomy_to_compare.split()[0].capitalize()
    in_species = args.taxonomy_to_compare.split()[1].lower().capitalize()
    genus = in_genus
    species = in_species
    total_tax = f"{genus} {species}    (selected manually)"

# Initialize variables for matching
found = False

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

            with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
                outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: {sample_gc_percent}")

            found = True

            break

# Handle unmatched cases
if not found:

    logging.INFO(f"No match found for '{genus} {species}'")

    with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: -1")

    with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")
    exit(0)

# Calculate and print ratio
if expected_length == "NA" or not expected_length:

    logging.INFO("No expected length was found to compare to")

    with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: NA\nIsolate_St.Devs: NA\nActual_length: {assembly_length}\nExpected_length: NA\nRatio: -1")

    with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: No Match Found\nSpecies_GC_Min: No Match Found\nSpecies_GC_Max: No Match Found\nSpecies_GC_Mean: No Match Found\nSpecies_GC_Count: No Match Found\nSample_GC_Percent: No Match Found")
    exit(0)

elif assembly_length == "NA" or not assembly_length:

    logging.INFO("No assembly length was found to compare with")

    with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_StDev: {stdev}\nIsolate_St.Devs: NA\nActual_length: NA\nExpected_length: {expected_length}\nRatio: -2")

    with open(f"{args.sample_name}_GC_content_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_GC_StDev: {gc_stdev}\nSpecies_GC_Min: {gc_min}\nSpecies_GC_Max: {gc_max}\nSpecies_GC_Mean: {gc_mean}\nSpecies_GC_Count: {gc_count}\nSample_GC_Percent: NA")
    exit(0)

ratio = float(assembly_length) / float(expected_length)

logging.INFO(f"Actual - {assembly_length}")
logging.INFO(f"Expected - {expected_length}")
logging.INFO(f"Ratio - {ratio:.4f}")

with open(f"{args.sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
    outfile.write(f"Tax: {total_tax}\nNCBI_TAXID: {taxid}\nSpecies_St.Dev: {stdev}\nIsolate_St.Devs: {stdevs}\nActual_length: {assembly_length}\nExpected_length: {expected_length}\nRatio: {ratio:.4f}")
