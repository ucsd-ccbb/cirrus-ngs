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
num_threads=${10}
chromosome=${11}

#logging
mkdir -p $log_dir
log_file=$log_dir/'split.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software_dir=/shared/workspace/software
samtools=$software_dir/samtools/samtools-1.1/samtools
sambamba=$software_dir/sambamba/0.4.7/bin/sambamba
mkdir -p $workspace

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ] || [ ! -f $workspace/$fastq_end1$file_suffix.bai ]
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
    aws s3 cp $input_address/$fastq_end1$download_suffix.bai $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix
fi
##END_DOWNLOAD##


##SPLIT##
$samtools view -b $workspace/$fastq_end1$file_suffix chr$chromosome > \
    $workspace/$fastq_end1.$chromosome.bam

$sambamba index -t $num_threads $workspace/$fastq_end1.$chromosome.bam \
    $workspace/$fastq_end1.$chromosome.bam.bai
##END_SPLIT##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.$chromosome.bam*" --recursive
##END_UPLOAD##
