#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
file_name=$4
is_single_end=$5

workspace=/shared/workspace
# path to the fastqc app
fastqc=/shared/workspace/software/FastQC/fastqc
# directory for output
miRNA = /shared/workspace/data_archive/SmallRNASeq/SmallRNASeq


mkdir -p $miRNA/$file_name

# redirecting all output to a file
exec 1>>$miRNA/$project_name'.log'
exec 2>>$miRNA/$project_name'.log'

if [ "$is_single_end" == "YES" ]
then
echo "running fastqc..."
# $file_name/$file_name$file_suffix -- this is the fastq file
# $fastqc $workspace/$file_name/$file_name$file_suffix -o $workspace/$file_name/
$fastqc $miRNA/$file_name$file_suffix -o $miRNA/$file_name/
else
echo "running fastqc..."
$fastqc $miRNA/$file_name/$file_name"_R1"$file_suffix -o $miRNA/$file_name/
$fastqc $miRNA/$file_name/$file_name"_R2"$file_suffix -o $miRNA/$file_name/
fi

