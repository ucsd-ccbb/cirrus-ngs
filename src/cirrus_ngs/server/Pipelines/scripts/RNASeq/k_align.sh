#!/bin/bash
# This shell script uses kallisto alignment.

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
num_threads=${10}

#logging
mkdir -p $log_dir/$fastq_end1
log_file=$log_dir/$fastq_end1/'k_align.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

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


##ALIGNMENT##
if [ "$fastq_end2" == "NULL" ]
then
    $kallisto/kallisto quant -i $kallisto/kallisto_index -o $workspace/ -t $num_threads --single -l 50 -s 20 \
    $workspace/$fastq_end1$file_suffix
else
    $kallisto/kallisto quant -i $kallisto/kallisto_index -o $workspace/ -t $num_threads -l 50 \
    $workspace/$fastq_end1$file_suffix $workspace/$fastq_end2$file_suffix
fi
##END_ALIGNMENT##

## Rename files ##
mv $workspace/abundance.h5 $workspace/$fastq_end1.abundance.h5
mv $workspace/abundance.tsv $workspace/$fastq_end1.abundance.tsv
mv $workspace/run_info.json $workspace/$fastq_end1.run_info.json

##UPLOAD##
aws s3 cp $workspace $output_address/ --exclude "*" --include "*.json*" --include "*.h5*" --include "*.tsv*" --recursive
##END_UPLOAD##
