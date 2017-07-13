#!/bin/bash

input_address=$1
sample_file=$2
sample_dir=$3
log_dir=$4

sample_dir="$sample_dir/${$sample_file/\..*}"

exec 1>>$log_dir/download.log
exec 2>>$log_dir/download.log

echo "Downloading $sample_file ..."

if [ ! -f $sample_dir/$sample_file ]; then
    aws s3 cp $input_address/$sample_file $sample_dir
    if [[ $sample_dir/$sample_file == *.gz ]]; then
        echo "Unzipping $sample_file ..."
        gunzip $sample_dir/$sample_file
        echo "$sample_file has been unzipped"
    fi
fi

echo "Finished downloading $sample_file"
