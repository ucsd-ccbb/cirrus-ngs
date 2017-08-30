#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'make_UCSC_file.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -d $workspace/tags_$fastq_end1 ]
then
    mkdir -p $workspace/tags_$fastq_end1

    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #always download forward reads
    aws s3 cp $input_address/tags_$fastq_end1 $workspace/tags_$fastq_end1 --recursive
fi
##END_DOWNLOAD##


##FASTQC##
check_exit_status "$make_UCSC_file $workspace/tags_$fastq_end1 -o auto" $JOB_NAME $status_file
##END_FASTQC##


##UPLOAD##
aws s3 sync $workspace/tags_$fastq_end1 $output_address/tags_$fastq_end1 
##END_UPLOAD##
