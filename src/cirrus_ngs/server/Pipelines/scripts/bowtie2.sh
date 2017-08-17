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
log_file=$log_dir/'bowtie2.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software=/shared/workspace/software
bowtie=$software/bowtie2-2.3.2-legacy
samtools=$software/samtools/samtools-1.1/samtools
# fa_file=$software/bowtie_index/hairpin_human/hairpin_human.fa   # fa file as the reference genome
basename=hairpin_human
cut=.cut
sam=.sam


echo $bowtie
echo $basename

mkdir -p $workspace


##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$cut$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$cut$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$cut$file_suffix".gz"
    fi

    #always download forward reads
    aws s3 cp $input_address/$project_name/$fastq_end1$download_suffix $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$project_name/$fastq_end2$download_suffix $workspace/
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##


##BOWTIE 2 ALIGNMENT##

# go to the directory of the index file
cd $software/bowtie2_index/

# indexing a reference genome (already built, moved under software/bowtie2_index)
# cd ~
# $bowtie/bowtie2-build $fa_file $basename

if [ "$fastq_end2" == "NULL" ]
then
    # single end
    $bowtie/bowtie2 -x $basename -U $workspace/$fastq_end1$cut$file_suffix \
    -S $workspace/$fastq_end1$sam
else
    # paired end: using the same output file name ($fastq_end1.sam)
    $bowtie/bowtie2 -x $basename -1 $workspace/$fastq_end1$cut$file_suffix \
    -2 $workspace/$fastq_end2/$cut$file_suffix -S $workspace/$fastq_end1$sam
fi

##END BOWTIE 2 ##

# TODO: Samtool stats on output (sam)
$samtools stats $workspace/$fastq_end1$sam

##UPLOAD##
aws s3 cp $workspace $output_address/$project_name --exclude "*" --include "*.sam*" --recursive

# upload the txt file from samtool stats
aws s3 cp $workspace $output_address/$project_name --exclude "*" --include "*.txt*" --recursive
##END_UPLOAD##

