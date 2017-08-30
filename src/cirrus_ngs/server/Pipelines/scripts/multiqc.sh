#!/bin/bash

project_name=$1
file_suffix=$2  # not used, for shell script convention
root_dir=$3
fastq_end1=$4    # always NULL, for shell script convention
fastq_end2=$5    # always NULL, for shell script convention
input_address=$6    # this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   # this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    # not used, for shell script convention
all_samples=${10}

mkdir -p $log_dir
log_file=$log_dir/'multiqc.log'
exec 1>>$log_file
exec 2>>$log_file

# set workspace to the multiqc directory so that multiqc would search for files
# in this directory
workspace=$root_dir/$project_name/multiqc

mkdir -p $workspace

# Download fastqc zip files and text files from alignment
for file in $all_samples; do
    aws s3 cp $input_address/$file/ $workspace/ --exclude "*" --include "*_fastqc.zip" --recursive
    aws s3 cp $input_address/$file/ $workspace/ --exclude "*" --include "*.txt" --recursive

done

# Multiqc on all
#       fastqc.zip files
#       .txt files from alignment

$multiqc -f $workspace -o $workspace

# Upload the html file to s3
aws s3 cp $workspace/multiqc_report.html $output_address/multiqc_report.html

