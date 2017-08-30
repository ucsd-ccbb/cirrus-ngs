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
log_file=$log_dir/$fastq_end1/'ht_align.log'
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

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/
    fi
fi
##END_DOWNLOAD##

# Star align
if [ "$fastq_end2" == "NULL" ]
then
    $STAR --genomeDir $STAR_genome --readFilesIn $workspace/$fastq_end1$file_suffix --outFileNamePrefix $workspace/$fastq_end1.

else
    $STAR --genomeDir $STAR_genome --readFilesIn $workspace/$fastq_end1$file_suffix $workspace/$fastq_end2$file_suffix \
    --outFileNamePrefix $workspace/$fastq_end1.
fi

if [ ! -f $workspace/$fastq_end1."Aligned.out.bam" ]; then
   $samtools view -Sb $workspace/$fastq_end1."Aligned.out.sam" > $workspace/$fastq_end1."Aligned.out.bam"
fi

if [ ! -f $workspace/$fastq_end1."Aligned.out.sorted.bam" ]; then
   $samtools sort -m 2G -@ 4 $workspace/$fastq_end1."Aligned.out.bam" $workspace/$fastq_end1."Aligned.out.sorted"
fi

if [ ! -f $workspace/$fastq_end1."Aligned.out.sorted.bam.bai" ]; then
   $samtools index $workspace/$fastq_end1."Aligned.out.sorted.bam"
fi
# End star align

# Upload
aws s3 cp $workspace $output_address/ --exclude "*" --include "*.Aligned*" --exclude "*.sam*" --recursive
##END_UPLOAD##
