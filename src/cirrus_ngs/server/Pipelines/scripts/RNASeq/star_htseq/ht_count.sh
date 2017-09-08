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

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'ht_count.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    check_exit_status "aws s3 cp $input_address/$fastq_end1$file_suffix $workspace/" $JOB_NAME $status_file
fi
##END_DOWNLOAD##

# Count reads #

if [ ! -f $workspace/$fastq_end1"_counts.txt" ]; then
   check_exit_status "$samtools view $workspace/$fastq_end1$file_suffix | sort -s -k 1,1 - \
   |  htseq-count - $STAR_ref_genes > $workspace/$fastq_end1"_counts.txt"" $JOB_NAME $status_file
fi

##UPLOAD##
check_exit_status "aws s3 cp $workspace $output_address/ --exclude "*" --include "*_counts.txt" --recursive" $JOB_NAME $status_file
##END_UPLOAD##
