#!/bin/bash
yes | sudo yum install perl-Env

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
num_threads=${11}   # number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'cal_expression.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace
trim=.trim

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

# RSEM calculate expression

if [ "$fastq_end2" == "NULL" ]
then
    $rsem --star --star-path $star_path -p $num_threads $workspace/$fastq_end1$trim$file_suffix \
    $STAR_index $workspace/$fastq_end1

else
    # paired end
    $rsem --star --star-path $star_path -p $num_threads --paired-end $workspace/$fastq_end1$trim$file_suffix \
    $workspace/$fastq_end2$trim$file_suffix $STAR_index $workspace/$fastq_end1
fi

# TODO: to be tested -- perform samtools stats, for multiqc purposes
check_exit_status "$samtools stats $workspace/$fastq_end1.transcript.bam > $workspace/$fastq_end1.txt" $JOB_NAME $status_file
if [ -f $workspace/$fastq_end1.txt ]
then
    echo "Finished samtools stats"
else
    echo "Failed samtools stats"
fi

##UPLOAD##
aws s3 cp $workspace $output_address/ --recursive
##END_UPLOAD##