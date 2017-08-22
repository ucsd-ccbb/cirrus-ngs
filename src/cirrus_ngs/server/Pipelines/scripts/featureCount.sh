#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped

#logging
mkdir -p $log_dir
log_file=$log_dir/'featureCount.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software=/shared/workspace/software
featureCounts=$software/featureCount/bin/featureCounts
annotation=$software/annotation_file/Homo_sapiens.GRCh37.75.gtf

sam=.sam

mkdir -p $workspace

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$sam ]
then
    #this is the suffix of the input from s3
    download_suffix=$sam

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

##Feature count##

if [ "$fastq_end2" == "NULL" ]
then
    # single end
    # TODO: right now using 5 threads
    $featureCounts -T 5 -a $annotation -t exon -g gene_id \
    -o $workspace/counts.txt $workspace/$fastq_end1$sam

else
    # paired end
    $featureCounts -p -a $annotation -t exon -g gene_id \
    -o $workspace/counts.txt $workspace/$fastq_end1$sam
fi

##END Feature count##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "*.txt*" --recursive
##END_UPLOAD##


