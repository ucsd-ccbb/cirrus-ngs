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
mkdir -p $log_dir
log_file=$log_dir/'fastqc.log'
curr_log_file=$log_dir/"fastqc_$fastq_end1.log"
exec 1>>$curr_log_file
exec 2>>$curr_log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

check_step_already_done $fastq_end1"_fastqc.html" $JOB_NAME $log_file $curr_log_file

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
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##


##FASTQC##
check_exit_status "$fastqc $workspace/$fastq_end1$file_suffix -o $workspace/" $fastq_end1"_fastqc.html" $JOB_NAME
mv $workspace/$fastq_end1$file_suffix"_fastqc.html" $workspace/$fastq_end1"_fastqc.html" 2>/dev/null
mv $workspace/$fastq_end1$file_suffix"_fastqc.zip" $workspace/$fastq_end1"_fastqc.zip" 2>/dev/null

if [ "$fastq_end2" != "NULL" ];
then
    check_exit_status "$fastqc $workspace/$fastq_end2$file_suffix -o $workspace/" $fastq_end1"_fastqc.html" $JOB_NAME
    mv $workspace/$fastq_end2$file_suffix"_fastqc.html" $workspace/$fastq_end2"_fastqc.html" 2>/dev/null
    mv $workspace/$fastq_end2$file_suffix"_fastqc.zip" $workspace/$fastq_end2"_fastqc.zip" 2>/dev/null
fi

##END_FASTQC##


##UPLOAD##
include_end1=$fastq_end1"_fastqc*"
include_end2=$fastq_end2"_fastqc*"
aws s3 cp $workspace $output_address --exclude "*" --include "$include_end1" --include "$include_end2" --recursive
##END_UPLOAD##

cat $curr_log_file >> $log_file
rm $curr_log_file
