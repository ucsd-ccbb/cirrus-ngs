#!/bin/sh

input_file=$1
output_file=$2

if [ ! -f $output_file ]; then
   /shared/workspace/software/FastQC/fastqc $input_file
fi
