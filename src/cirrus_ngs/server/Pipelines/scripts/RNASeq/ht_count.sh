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
mkdir -p $log_dir/$fastq_end1
log_file=$log_dir/$fastq_end1/'ht_count.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    aws s3 cp $input_address/$fastq_end1$file_suffix $workspace/
fi
##END_DOWNLOAD##

# Count reads #

if [ ! -f $workspace/$fastq_end1"_counts.txt" ]; then
   $samtools view $workspace/$fastq_end1$file_suffix | sort -s -k 1,1 - |  htseq-count - $ref_genes > $workspace/$fastq_end1"_counts.txt"
fi

##UPLOAD##
aws s3 cp $workspace $output_address/ --exclude "*" --include "*_counts.txt" --recursive
##END_UPLOAD##
