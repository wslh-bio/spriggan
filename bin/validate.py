#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse

#this gets us the root dir of the project
base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('spriggan_report_valid',help='Path to validated spriggan_report.csv')
parser.add_argument('spriggan_report_test',help='Path to test spriggan_report.csv')
args = parser.parse_args()

valid_results = pd.read_csv(os.path.abspath(args.spriggan_report_valid),sep=',',index_col="Sample").sort_index()
test_results = pd.read_csv(os.path.abspath(args.spriggan_report_test),sep=',',index_col="Sample").sort_index()

### Sort Columns By name
valid_results = valid_results.reindex(sorted(valid_results.columns),axis=1)
test_results = test_results.reindex(sorted(test_results.columns),axis=1)

### Validate Results
validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If no difference validation is successful
if validation.empty:
    print("Validation check Successful!")
    sys.exit()

### If assembly length differs by less than 1000 bp then remove from dataframe
if "Assembly Length (bp)" in validation.columns:
    for sample in validation["Assembly Length (bp)"].index.get_level_values('Sample').unique():
        valid_data = validation["Assembly Length (bp)"].loc[sample,"Valid Data"]
        test_data = validation["Assembly Length (bp)"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        if diff < 1000:
            test_results.loc[sample,"Assembly Length (bp)"] = valid_results.loc[sample,"Assembly Length (bp)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If contig number differs by less than 10 then set test equal to valid
if "Contigs (#)" in validation.columns:
    for sample in validation["Contigs (#)"].index.get_level_values('Sample').unique():
        valid_data = validation["Contigs (#)"].loc[sample,"Valid Data"]
        test_data = validation["Contigs (#)"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        if diff < 10:
            test_results.loc[sample,"Contigs (#)"] = valid_results.loc[sample,"Contigs (#)"]
            validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

if "AMR" in validation.columns:
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

### If no difference validation is successful
if validation.empty:
    print("Validation check Successful!")
    sys.exit()
else:
    print("Validation Failed")
    print(validation)
    sys.exit(1)