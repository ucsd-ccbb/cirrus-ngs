#!/bin/bash

input_dir=$1
sample_dir=$2
file1_name=$3
file2_name=$4
log_dir=$5


exec 1>>$log_dir/fastqc.log
exec 2>>$log_dir/fastqc.log

sample_dir="$sample_dir/fastqc"
file1_name=${file1_name//.gz/}
file2_name=${file2_name//.gz/}

output1="${file1_name}_fastqc.zip"
output2="${file2_name}_fastqc.zip"

mkdir -p $sample_dir

if [ ! -f $sample_dir/$output1 ]; then
    echo "Performing FastQC analysis on $file1_name ..."

    /shared/workspace/software/FastQC/fastqc --noextract -o $sample_dir \
        $input_dir/$file1_name

    echo "Finished FastQC analysis on $file1_name"
else
    echo "FastQC analysis has already been performed on $file1_name"
fi
echo

if [ "$file2_name" != "NULL" ]; then
    if [ ! -f $sample_dir/$output2 ]; then
        echo "Performing FastQC analysis on $file2_name ..."

        /shared/workspace/software/FastQC/fastqc --noextract -o $sample_dir \
            $input_dir/$file2_name

        echo "Finished FastQC analysis on $file2_name"
    else
        echo "FastQC analysis has already been performed on $file2_name"
    fi
fi

##DEBUG##
echo
echo "dir is $sample_dir"
echo `ls $sample_dir`
##ENDDEBUG##
echo
