#!/bin/bash

s3_input_files_address=$1
sample_file=$2
sample_dir=$3

if [ ! -f $3/$2 ]; then
   aws s3 cp $1/$2 $3
fi
