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

sam=.sam

#logging
mkdir -p $log_dir
log_file=$log_dir/'counting.log'
exec 1>>$log_file
exec 2>>$log_file

workspace=$root_dir/$project_name/$workflow
mkdir -p $workspace

# Download files from s3
for file in $all_samples; do
    if [ ! -f $workspace/$file$sam ]
    then
        aws s3 cp $input_address/$file/$file$sam $workspace/
    fi
done

# Call counter.count
check_exit_status "/shared/workspace/software/anaconda3/bin/python /shared/workspace/Pipelines/util/counter.py \
$hairpin_human_fa $workspace/counts.out $workspace/rates.out "$all_samples" $workspace" $JOB_NAME $status_file

# Upload the output file to s3
aws s3 cp $workspace $output_address --exclude "*" --include "*.out*" --recursive
