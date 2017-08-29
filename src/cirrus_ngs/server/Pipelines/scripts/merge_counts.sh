#!/bin/bash

project_name=$1
file_suffix=$2
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9
all_samples=${10}    # a space separated string containing all the sample names
num_threads=${11}
workflow=${12}

#logging
mkdir -p $log_dir
log_file=$log_dir/'merge.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name
mkdir -p $workspace

# Download files from s3
for file in $all_samples; do
    if [ ! -f $workspace/$file$file_suffix ]
    then
        aws s3 cp $input_address/$file/$file$file_suffix $workspace/
    fi
done

# Call the merge count file
python /shared/workspace/Pipelines/util/MergeCountFile.py $workflow "$all_samples" $workspace

# Upload the output file
aws s3 cp $workspace $output_address/ --exclude "*" --include "*all_gene_counts.txt*" --recursive