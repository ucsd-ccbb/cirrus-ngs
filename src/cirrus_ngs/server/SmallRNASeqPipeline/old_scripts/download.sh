#!/bin/bash

input_address=$1
sample_dir=$2
file1_name=$3
file2_name=$4
logs_dir=$5

sample_dir="$sample_dir/${file1_name//.*/}"

exec 1>>$logs_dir/download.log
exec 2>>$logs_dir/download.log

echo "Downloading $file1_name ..."

if [ ! -f $sample_dir/$file1_name ]; then
aws s3 cp $input_address/$file1_name $sample_dir
if [[ $sample_dir/$file1_name == *.gz ]]; then
echo "Unzipping $file1_name ..."
gunzip $sample_dir/$file1_name
echo "$file1_name has been unzipped"
fi
fi

echo "Finished downloading $file1_name"

if [ "$file2_name" != "NULL" ]; then
echo "Downloading $file2_name ..."
aws s3 cp $input_address/$file2_name $sample_dir
if [[ $sample_dir/$file2_name == *gz ]]; then
echo "Unzipping $file2_name ..."
gunzip $sample_dir/$file2_name
echo "$file2_name has been unzipped"
fi

echo "Finished downloading $file2_name"
fi

echo `ls $sample_dir`

