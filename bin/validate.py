#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse
import logging
import sys

#Setting up logging structure
logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

#this gets us the root dir of the project
base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.', 
                                 epilog='Example usage: python validate.py $VALID_REPORT $TEST_REPORT $VALID_STDEV_FILE')
parser.add_argument('spriggan_report_valid',
                    help='Path to validated spriggan_report.csv')
parser.add_argument('spriggan_report_test',
                    help='Path to test spriggan_report.csv')
parser.add_argument('valid_stdev_file',
                    help='Path to validate calculate folder')

args = parser.parse_args()

valid_results = pd.read_csv(os.path.abspath(args.spriggan_report_valid),sep=',',index_col="Sample").sort_index()
test_results = pd.read_csv(os.path.abspath(args.spriggan_report_test),sep=',',index_col="Sample").sort_index()
stdev = pd.read_csv(os.path.abspath(args.valid_stdev_file),sep=',',index_col="Sample").sort_index()

### Sort Columns By name
valid_results = valid_results.reindex(sorted(valid_results.columns),axis=1)
test_results = test_results.reindex(sorted(test_results.columns),axis=1)
stdev = stdev.reindex(sorted(stdev.columns),axis=1)

### Validate Results
validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If no difference validation is successful
if validation.empty:
    logging.info("Validation check Successful!")
    sys.exit()

### If assembly length differs by less than 1000 bp then remove from dataframe
if "Assembly Length (bp)" in validation.columns:

    logging.debug("Processing assembly length.")

    for sample in validation["Assembly Length (bp)"].index.get_level_values('Sample').unique():
        valid_data = validation["Assembly Length (bp)"].loc[sample,"Valid Data"]
        test_data = validation["Assembly Length (bp)"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)

        logging.debug("Processing difference of assembly length")
        if diff < 1000:
            test_results.loc[sample,"Assembly Length (bp)"] = valid_results.loc[sample,"Assembly Length (bp)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If contig number differs by less than 10 then set test equal to valid
if "Contigs (#)" in validation.columns:
    
    logging.debug("Processing Contig numbers.")
    
    for sample in validation["Contigs (#)"].index.get_level_values('Sample').unique():
        valid_data = validation["Contigs (#)"].loc[sample,"Valid Data"]
        test_data = validation["Contigs (#)"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        
        logging.debug("Processing difference of contig numbers")

        if diff < 10:
            test_results.loc[sample,"Contigs (#)"] = valid_results.loc[sample,"Contigs (#)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

if "AMR" in validation.columns:
    
    logging.debug("Processing AMR genes.")

    for sample in validation["AMR"].index.get_level_values('Sample').unique():
        passing = True
        valid_data_mech = validation["AMR"].loc[sample,"Valid Data"].split(';')
        test_data_mech = validation["AMR"].loc[sample,"Test Data"].split(';')
        valid_data_coverage = valid_results.loc[sample,"AMR Coverage"].split(';')
        test_data_coverage = test_results.loc[sample,"AMR Coverage"].split(';')
        valid_data_identity = valid_results.loc[sample,"AMR Identity"].split(';')
        test_data_identity = test_results.loc[sample,"AMR Identity"].split(';')

        valid_data = sorted(list(zip(valid_data_mech,valid_data_coverage,valid_data_identity)), key=lambda x: x[0])
        test_data = sorted(list(zip(test_data_mech,test_data_coverage,test_data_identity)), key=lambda x: x[0])

        for v,t in zip(valid_data,test_data):
            if v != t:
                passing = False
        
        if passing:
            test_results.loc[sample,"AMR"] = valid_results.loc[sample,"AMR"]
            test_results.loc[sample,"AMR Coverage"] = valid_results.loc[sample,"AMR Coverage"]
            test_results.loc[sample,"AMR Identity"] = valid_results.loc[sample,"AMR Identity"]
    validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

if "Selected AMR Genes" in validation.columns:
    
    logging.debug("Processing selected AMR genes.")
    
    for sample in validation["Selected AMR Genes"].index.get_level_values('Sample').unique():
        passing = True
        valid_data_mech = validation["Selected AMR Genes"].loc[sample,"Valid Data"].split(';')
        test_data_mech = validation["Selected AMR Genes"].loc[sample,"Test Data"].split(';')
        valid_data_coverage = valid_results.loc[sample,"Selected AMR Genes Coverage"].split(';')
        test_data_coverage = test_results.loc[sample,"Selected AMR Genes Coverage"].split(';')
        valid_data_identity = valid_results.loc[sample,"Selected AMR Genes Identity"].split(';')
        test_data_identity = test_results.loc[sample,"Selected AMR Genes Identity"].split(';')

        valid_data = sorted(list(zip(valid_data_mech,valid_data_coverage,valid_data_identity)), key=lambda x: x[0])
        test_data = sorted(list(zip(test_data_mech,test_data_coverage,test_data_identity)), key=lambda x: x[0])

        for v,t in zip(valid_data,test_data):
            if v != t:
                passing = False
        
        if passing:
            test_results.loc[sample,"Selected AMR Genes"] = valid_results.loc[sample,"Selected AMR Genes"]
            test_results.loc[sample,"Selected AMR Genes Coverage"] = valid_results.loc[sample,"Selected AMR Genes Coverage"]
            test_results.loc[sample,"Selected AMR Genes Identity"] = valid_results.loc[sample,"Selected AMR Genes Identity"]
    validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### Checks to see if MLST Scheme is in different order than validated script
if "MLST Scheme" in validation.columns:

    for sample in validation["MLST Scheme"].index.get_level_values('Sample').unique():

        logging.debug("Setting up comparison")
        valid_data = validation["MLST Scheme"].loc[sample,"Valid Data"]
        test_data = validation["MLST Scheme"].loc[sample,"Test Data"]

        try:
            logging.debug("Sorting values if ; is present")
            if ";" in test_data and ";" in valid_data:
                valid_data = sorted(valid_data.split(";"))
                test_data = sorted(test_data.split(";"))
        except TypeError:
            logging.debug("If no ; is present and test data = valid data, pass.")
            if test_data == valid_data:
                pass

        logging.debug("Updating if passing")
        if valid_data == test_data:
            passing = True
        else:
            passing = False

        logging.debug("If passing adding to Valid results")
        if passing:
            test_results.loc[sample,"MLST Scheme"] = valid_results.loc[sample, "MLST Scheme"]

    logging.debug("Adding MLST to valid results")
    validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### Checks if genome length ratio is within 1 standard deviation from the mean
if "Genome Length Ratio (Actual/Expected)" in validation.columns:
    logging.debug("Processing genome length ratios.")

    for sample in validation["Genome Length Ratio (Actual/Expected)"].index.get_level_values('Sample').unique():

        valid_data = validation["Genome Length Ratio (Actual/Expected)"].loc[sample,"Valid Data"]
        test_data = validation["Genome Length Ratio (Actual/Expected)"].loc[sample,"Test Data"]
        assembly_stdev = stdev.loc[sample, 'assembly_stdev']

        logging.debug("Calculate lower and higher bounds of Genome Length Ratio.")

        lower = valid_data - assembly_stdev
        higher = valid_data + assembly_stdev
        
        logging.debug("Check if test_data is within the bounds for genome length ratio.")

        if lower <= test_data <= higher:
            test_results.loc[sample,"Genome Length Ratio (Actual/Expected)"] = valid_results.loc[sample,"Genome Length Ratio (Actual/Expected)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### Checks if sample GC content is within 1 standard deviation from the species GC mean
if "Sample GC Content (%)" in validation.columns:
    logging.debug("Process sample GC content.")

    for sample in validation["Sample GC Content (%)"].index.get_level_values('Sample').unique():

        valid_sample_data = validation["Sample GC Content (%)"].loc[sample,"Valid Data"]
        test_sample_data = validation["Sample GC Content (%)"].loc[sample,"Test Data"]
        gc_stdev = stdev.loc[sample, 'species_gc_stdev']
        gc_mean = valid_results.loc[sample, 'Species GC Content (Mean)']

        logging.debug("Calculate lower and higher bounds.")

        lower = gc_mean - gc_stdev
        higher = gc_mean + gc_stdev 

        logging.debug("Check if test_data is within the bounds.")

        if lower <= test_sample_data <= higher:

            test_results.loc[sample,"Sample GC Content (%)"] = valid_results.loc[sample,"Sample GC Content (%)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If no difference validation is successful
if validation.empty:
    logging.info("Validation check Successful!")
    sys.exit()
else:
    logging.info("Validation Failed")
    logging.info(validation)
    sys.exit(1)