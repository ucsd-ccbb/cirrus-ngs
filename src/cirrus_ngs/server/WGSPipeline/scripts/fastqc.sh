#!/bin/bash

input=$1
output=$2
log_dir=$3

exec 1>>$log_dir/fastqc.log
exec 2>>$log_dir/fastqc.log

if [ ! -d $output ]; then
    mkdir $output
fi

/shared/workspace/software/FastQC/fastqc --noextract -o $output $input
