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
platform_technoogy=${11}

#logging
mkdir -p $log_dir
log_file=$log_dir/'bwa.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software_dir=/shared/workspace/software
genome=$software_dir/genomes/Hsapiens/bwa/ucsc.hg19.fasta
bwa=$software_dir/bwa/bwa-0.7.12/bwa
samblaster=$software_dir/samblaster/samblaster
samtools=$software_dir/samtools/samtools-1.1/samtools
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


##ALIGN##
if [ "fastq_end2" == "NULL" ]
then
    $bwa mem -M -t $num_threads -R "@RG\tID:1\tPL:ILLUMINA\tPU:tempID\tSM:$fastq_end1" -v 1 \
        $genome $workspace/$fastq_end1$file_suffix | $samblaster | \
        $samtools view -Sb - > $workspace/$fastq_end1.bam
else
    $bwa mem -M -t $num_threads -R "@RG\tID:1\tPL:ILLUMINA\tPU:tempID\tSM:$fastq_end1" -v 1 \
        $genome $workspace/$fastq_end1$file_suffix \
        $workspace/$fastq_end2$file_suffix | $samblaster | \
        $samtools view -Sb - > $workspace/$fastq_end1.bam
fi
##END_ALIGN##

##UPLOAD##
aws s3 cp $workspace $output_address/ --exclude "*" --include "$fastq_end1.bam" --recursive
##END_UPLOAD##
