#!/bin/bash

project_name=$1
file_suffix=$2
file_path=$3
file1_name=$4
file2_name=$5
input_address=$6

workspace=$file_path/$project_name

if [ ! -d $workspace/$file1_name ]; then
   mkdir -p $workspace/$file1_name
fi

# redirecting all output to a file
exec 1>>$workspace/'download.log'
exec 2>>$workspace/'download.log'

if [ $file2_name == "NULL" ]; then
    echo "Downloading $file1_name ..."
    aws s3 cp $input_address/$file1_name$file_suffix".gz" $workspace/$file1_name/
    gunzip $workspace/$file1_name/$file1_name$file_suffix".gz"

    echo "Finished downloading $file1_name"
else
    echo "Downloading $file1_name and $file2_name..."
    aws s3 cp $input_address/$file1_name$file_suffix".gz" $workspace/$file1_name/
    gunzip $workspace/$file1_name/$file1_name$file_suffix".gz"
    aws s3 cp $input_address/$file2_name$file_suffix".gz" $workspace/$file1_name/
    gunzip $workspace/$file1_name/$file2_name$file_suffix".gz"

    echo "Finished downloading $file1_name and $file2_name"
fi

