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
log_file=$log_dir/'fastqc.log'
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
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #always download forward reads
    check_exit_status "aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/" $JOB_NAME $status_file
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        check_exit_status "aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/" $JOB_NAME $status_file
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##


##FASTQC##
check_exit_status "$fastqc $workspace/$fastq_end1$file_suffix -o $workspace/" $JOB_NAME $status_file
mv $workspace/$fastq_end1$file_suffix"_fastqc.html" $workspace/$fastq_end1"_fastqc.html" 2>/dev/null
mv $workspace/$fastq_end1$file_suffix"_fastqc.zip" $workspace/$fastq_end1"_fastqc.zip" 2>/dev/null

if [ "$fastq_end2" != "NULL" ];
then
    check_exit_status "$fastqc $workspace/$fastq_end2$file_suffix -o $workspace/" $JOB_NAME $status_file
    mv $workspace/$fastq_end2$file_suffix"_fastqc.html" $workspace/$fastq_end2"_fastqc.html" 2>/dev/null
    mv $workspace/$fastq_end2$file_suffix"_fastqc.zip" $workspace/$fastq_end2"_fastqc.zip" 2>/dev/null
fi
##END_FASTQC##


##UPLOAD##
include_end1=$fastq_end1"_fastqc*"
include_end2=$fastq_end2"_fastqc*"
check_exit_status "aws s3 cp $workspace $output_address --exclude "*" --include "$include_end1" --include "$include_end2" --recursive" $JOB_NAME $status_file
##END_UPLOAD##
