#!/bin/bash

project_name=$1
file_suffix=$2
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9

#logging
mkdir -p $log_dir/$fastq_end1
log_file=$log_dir/$fastq_end1/'k_count.log'
exec 1>>$log_file
exec 2>>$log_file

echo $output_address

# sample=$1
# output_file=$2

# Call the perl file
perl $kallisto_counts $output_address



