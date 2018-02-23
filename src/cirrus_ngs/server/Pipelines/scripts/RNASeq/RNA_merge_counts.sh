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
all_samples=${11}    # a space separated string containing all the sample names

#logging
mkdir -p $log_dir
log_file=$log_dir/'merge_counts.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

# Download files from s3
for file in $all_samples; do
    check_exit_status "aws s3 cp $input_address/$file/$file$file_suffix $workspace/" $JOB_NAME $status_file
done

# Call the merge count file
check_exit_status "python /shared/workspace/Pipelines/util/RNA_MergeCount.py $workflow $workspace $all_samples" $JOB_NAME $status_file

# Upload the output file
aws s3 cp $workspace $output_address/ --exclude "*" --include "all_gene_counts.txt" \
    --recursive
