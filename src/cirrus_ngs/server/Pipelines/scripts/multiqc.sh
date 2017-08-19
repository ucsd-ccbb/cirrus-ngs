#!/bin/bash

project_name=$1
root_dir=$2            # /scratch
input_address=$3    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$4   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$5


mkdir -p $log_dir
log_file=$log_dir/'multiqc.log'
exec 1>>$log_file
exec 2>>$log_file

# set workspace to the multiqc directory so that all the files would be
# the input for multiqc tool
workspace=$root_dir/$project_name/multiqc
software=/shared/workspace/software
multiqc=$software/anaconda3/bin/multiqc

mkdir -p $workspace

# Download fastqc zip files
aws s3 cp $input_address $workspace --exclude "*" --include "*_fastqc.zip" --recursive
# Download text files from the output of cutadapt and alignment (text generated during the run of that tool)
aws s3 cp $input_address $workspace --exclude "*" --include "*.txt" --recursive

# Multiqc on all
#       .txt files;
#       fastqc.zip files (and more)
$multiqc -f $workspace --cl_config "sp: {cutadapt: [fn: 'cut.txt', fn: 'cut_3.txt']}" -o $workspace

# Upload the html file to s3
aws s3 cp $workspace/multiqc_report.html $output_address/multiqc_report.html

