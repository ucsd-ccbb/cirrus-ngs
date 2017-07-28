#!/bin/bash

input_address=$1
sample_dir=$2
file1_name=$3
file2_name=$4
log_dir=$5

sample_dir="$sample_dir/${file1_name//.*/}"

exec 1>>$log_dir/download.log
exec 2>>$log_dir/download.log

if [ ! -f `sed s/"\.gz"// <<< "$sample_dir/$file1_name"` ]; then
    echo "Downloading $file1_name ..."
    aws s3 cp $input_address/$file1_name $sample_dir
    if [[ $sample_dir/$file1_name == *.gz ]]; then
        echo "Unzipping $file1_name ..."
        gunzip $sample_dir/$file1_name
        echo "$file1_name has been unzipped"
    fi
    echo "Finished downloading $file1_name"
else
    echo "$file1_name has already been downloaded"
fi
echo


if [ "$file2_name" != "NULL" ]; then
    if [ ! -f `sed s/"\.gz"// <<< "$sample_dir/$file2_name"` ]; then
        echo "Downloading $file2_name ..."
        aws s3 cp $input_address/$file2_name $sample_dir
        if [[ $sample_dir/$file2_name == *.gz ]]; then
            echo "Unzipping $file2_name ..."
            gunzip $sample_dir/$file2_name
            echo "$file2_name has been unzipped"
        fi
        echo "Finished downloading $file2_name"
    else
        echo "$file2_name has already been downloaded"
    fi
fi

##DEBUG##
echo
echo "dir is $sample_dir"
echo `ls $sample_dir`
##ENDDEBUG##
echo
