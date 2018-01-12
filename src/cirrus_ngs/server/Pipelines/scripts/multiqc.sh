#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
fastq_end1=$5
fastq_end2=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
all_samples=${11}

mkdir -p $log_dir
log_file=$log_dir/'multiqc.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

# set workspace to the multiqc directory so that multiqc would search for files
# in this directory
workspace=$root_dir/$project_name/multiqc
rm -r $workspace &>/dev/null
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

# Download fastqc zip files and text files generated with the alignment
for file in $all_samples
do
    aws s3 cp $input_address/$file/ $workspace/ --exclude "*" --include "*_fastqc.zip" --recursive
    aws s3 cp $input_address/$file/ $workspace/ --exclude "*" --include "$file.txt" --recursive
done

# Multiqc on all
#       .fastqc.zip files
#       .txt files from alignment

check_exit_status "$multiqc -f $workspace -o $workspace" $JOB_NAME $status_file

# Upload the html file to s3
aws s3 cp $workspace/multiqc_report.html $output_address/multiqc_report.html
