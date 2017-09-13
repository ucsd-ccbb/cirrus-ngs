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
num_threads=${11}
min_len=${12}       # drop the read if it is below this minimum length

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'cutadapt.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace
adapter=TGGAATTCTCGGGTGCCAAGG
trim=.trim
cut=.cut

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file


##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$trim$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$trim$file_suffix

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

##CUT_ADAPT##
# cut 3' end
check_exit_status "$cutadapt -a $adapter -o $workspace/$fastq_end1$cut$file_suffix \
$workspace/$fastq_end1$trim$file_suffix -m $min_len" $JOB_NAME $status_file

if [ "$fastq_end2" != "NULL" ];
then
    # cut 3' end
    check_exit_status "$cutadapt -a $adapter -o $workspace/$fastq_end2$cut$file_suffix \
    $workspace/$fastq_end2$trim$file_suffix -m $min_len" $JOB_NAME $status_file
fi
##END_CUT_ADAPT##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "*.cut.fastq*" --recursive
##END_UPLOAD##
