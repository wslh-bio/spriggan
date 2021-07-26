#!/bin/bash

# AR Summary Check
#ar_test=$( grep -vFxf ./autocheck/ar_summary_std.tsv $1 )
#if [ ! -z "$ar_test" ]
#then
#    echo "Failed AR Check:"
#    echo $ar_test
#    exit 1
#fi

# Coverage Check
#temp_cov_file=$( mktemp )
#awk 'NR > 1 { $2 = sprintf("%s%d", "\t", $2 + 0.5); $3 = sprintf("%s%d", "\t", $3 + 0.5) }; 1' $2 > $temp_cov_file
#cov_test=$( grep -vFxf ./autocheck/coverage_stats_std.tsv $temp_cov_file )
#if [ ! -z "$cov_test" ]
#then
#    echo "Failed Coverage Check:"
#    echo $cov_test
#    exit 1
#fi
#rm $temp_cov_file

# Kraken Check
krak_test=$( grep -vFxf ./autocheck/kraken_results_std.txt $3 )
if [ ! -z "$krak_test" ]
then
    echo "Failed Kraken Check:"
    echo $krak_test
    exit 1
fi

# MLST Check
mlst_test=$( grep -vFxf ./autocheck/mlst_formatted_std.tsv $4 )
if [ ! -z "$mlst_test" ]
then
    echo "Failed MLST Check:"
    echo $mlst_test
    exit 1
fi

# QUAST Check
#quast_test=$( grep -vFxf ./autocheck/quast_results_std.tsv $5 )
#if [ ! -z "$quast_test" ]
#then
#    echo "Failed QUAST Check:"
#    echo $quast_test
#    exit 1
#fi
