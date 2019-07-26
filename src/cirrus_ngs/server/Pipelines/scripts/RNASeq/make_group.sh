#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
fastq_end1=$5
fastq_end2=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped

#logging
mkdir -p $log_dir
log_file=$log_dir/'make_group.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

# Make a new group file
touch $workspace/group.txt

# Call the group file maker
$python /shared/workspace/cirrus-ngs/src/cirrus_ngs/server/Pipelines/util/GroupFileMaker.py $workspace/group.txt \
    /shared/workspace/cirrus-ngs/src/cirrus_ngs/server/Pipelines/yaml_files/RNASeq/$workflow/$project_name.yaml $output_address

# View the group file
head $workspace/group.txt

# Upload the group file
aws s3 cp $workspace $output_address/ --exclude "*" --include "*group.txt*" --recursive --quiet
