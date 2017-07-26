#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
file1_name=$4
file2_name=$5
genome=$6

software_dir=/shared/workspace/software
workspace=$file_path/$project_name
fastqc=$software_dir/FastQC/fastqc

if [ ! -d $workspace/$file1_name ]; then
   mkdir -p $workspace/$file1_name
fi

# redirecting all output to a file
exec 1>>$workspace/'fastqc.log'
exec 2>>$workspace/'fastqc.log'

if [ "$file2_name" == "NULL" ];
then
   echo “Running fastqc on $file1_name...”
   $fastqc $workspace/$file1_name/$file1_name$file_suffix -o $workspace/$file1_name/

   echo "Finished fastqc on $file1_name"
else
   echo “Running fastqc on $file1_name and $file2_name...”
   $fastqc $workspace/$file1_name/$file1_name$file_suffix -o $workspace/$file1_name/
   $fastqc $workspace/$file1_name/$file2_name$file_suffix -o $workspace/$file1_name/

   echo "Finished fastqc on $file1_name and $file2_name"
fi

