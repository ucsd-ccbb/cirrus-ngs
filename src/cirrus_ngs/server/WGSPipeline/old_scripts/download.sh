#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6
is_zipped=$7

workspace=$file_path/$project_name

mkdir -p $workspace/$fastq_end1
mkdir -p ~/$project_name/

# redirecting all output to a file
exec 1>>~/$project_name/download.log
exec 2>>~/$project_name/download.log

if [ "$is_zipped" == "True" ];
then
    file_suffix=$file_suffix".gz"
fi

echo "Downloading $fastq_end1 ..."

if [ $fastq_end2 == "NULL" ]; 
then
    aws s3 cp $input_address/$fastq_end1$file_suffix $workspace/$fastq_end1/
    gunzip $workspace/$fastq_end1/$fastq_end1$file_suffix
else
    aws s3 cp $input_address/$fastq_end1$file_suffix $workspace/$fastq_end1/
    gunzip $workspace/$fastq_end1/$fastq_end1$file_suffix
    aws s3 cp $input_address/$fastq_end2$file_suffix $workspace/$fastq_end1/
    gunzip $workspace/$fastq_end1/$fastq_end2$file_suffix
fi

echo 
echo `ls $workspace/$fastq_end1`
