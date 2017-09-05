#!/bin/bash

project_name=$1
file_suffix=$2  # extension of input file
root_dir=$3
fastq_end1=$4    # always NULL, for shell script convention
fastq_end2=$5    # always NULL, for shell script convention
input_address=$6    # this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   # this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    # always "False" in this case
all_samples=${10}

fa_file=/shared/workspace/software/bowtie_index/hairpin_human/hairpin_human.fa
sam=.sam

mkdir -p $log_dir
log_file=$log_dir/'counting.log'
exec 1>>$log_file
exec 2>>$log_file

workspace=$root_dir/$project_name

mkdir -p $workspace

# Download files from s3
for file in $all_samples; do
    if [ ! -f $workspace/$file$sam ]
    then
        aws s3 cp $input_address/$file/$file$sam $workspace/
    fi
done

# Call counter.count
/shared/workspace/software/anaconda3/bin/python /shared/workspace/Pipelines/util/counter.py \
$fa_file $workspace/counts.out $workspace/rates.out "$all_samples" $workspace

# Upload the output file to s3
aws s3 cp $workspace $output_address --exclude "*" --include "*.out*" --recursive
