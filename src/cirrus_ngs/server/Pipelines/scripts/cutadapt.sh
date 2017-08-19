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
log_file=$log_dir/'cutadapt.log'
exec 1>>$log_file
exec 2>>$log_file

echo $project_name
echo $file_suffix
echo $root_dir
echo $fastq_end1
echo $fastq_end2
echo $input_address
echo $output_address
echo $log_dir
echo $is_zipped

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software=/shared/workspace/software
cutadapt=$software/anaconda3/bin/cutadapt
adapter=TGGAATTCTCGGGTGCCAAGG
trim=.trim
cut=.cut
threeprime=_3

mkdir -p $workspace


##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$trim$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$trim$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$trim$file_suffix.gz
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


##CUT_ADAPT##
# TODO: for paired end too
# cut 3' end
$cutadapt -a $adapter -o $workspace/$fastq_end1$cut$threeprime$file_suffix \
$workspace/$fastq_end1$trim$file_suffix -m 2 --info-file $workspace/$fastq_end1$cut$threeprime.txt
# cut anchored 5' end
$cutadapt -g ^$adapter -o $workspace/$fastq_end1$cut$file_suffix \
$workspace/$fastq_end1$cut$threeprime$file_suffix -m 2 --info-file $workspace/$fastq_end1$cut.txt

if [ "$fastq_end2" != "NULL" ];
then
    # cut 3' end
    $cutadapt -a $adapter -o $workspace/$fastq_end2$cut$threeprime$file_suffix \
    $workspace/$fastq_end2$trim$file_suffix -m 2
    # cut anchored 5' end
    $cutadapt -g ^$adapter -o $workspace/$fastq_end2$cut$file_suffix \
    $workspace/$fastq_end2$cut$threeprime$file_suffix -m 2
fi
##END_CUT_ADAPT##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "*.cut.fastq*" --recursive

# upload the txt files
aws s3 cp $workspace $output_address --exclude "*" --include "*.txt*" --recursive
##END_UPLOAD##

