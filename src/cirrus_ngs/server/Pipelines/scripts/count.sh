#!/bin/bash

project_name=$1
root_dir=$2            # /scratch
input_address=$3    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$4   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$5
fa_file=$6
counts_out=$7
mapping_rates_out=$8
array=( $@ )
len=${#array[@]}
samfiles_list=${array[@]:8:$len}  # all samfiles as one comma separated string

mkdir -p $log_dir
log_file=$log_dir/'counting.log'
exec 1>>$log_file
exec 2>>$log_file

echo $log_dir
echo $fa_file
echo $counts_out
echo $mapping_rates_out
echo $samfiles_list

workspace=$root_dir/$project_name

mkdir -p $workspace

# Download files from s3
for i in ${samfiles_list//,/ }; do
    echo "In shell: "$i
    if [ ! -f $workspace/$i ]
    then
        aws s3 cp $input_address/$project_name/$i $workspace
    fi
done

# Call counter.count
/shared/workspace/software/anaconda3/bin/python /shared/workspace/Pipelines/util/counter.py \
$fa_file $workspace/$counts_out $workspace/$mapping_rates_out $samfiles_list $workspace

# Upload the output file to s3
aws s3 cp $workspace $output_address/$project_name --exclude "*" --include "*.out*" --recursive
