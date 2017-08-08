#!/bin/bash

project_name=$1
root_dir=$2            # /scratch
input_address=$3    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$4   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$5
notebook_dir = $6   # passed in from python, path to the notebook
array=( $@ )
len=${#array[@]}
multiqc_file_list=${array[@]:6:$len}  # multiqc fastqc.zip files: unknown number of args

echo $project_name
echo $root_dir
echo $input_address
echo $output_address
echo $log_dir
echo $notebook_dir
echo $multiqc_file_list

mkdir -p $log_dir
log_file=$log_dir/'multiqc.log'
exec 1>>$log_file
exec 2>>$log_file

workspace=$root_dir/$project_name
software=/shared/workspace/software
multiqc=$software/anaconda3/bin/multiqc

mkdir -p $workspace

# Download fastqc.zip files from s3
for i in "$@"; do
    if [ ! -f $workspace/$i ]
    then
        aws s3 cp $input_address/$i $workspace/
    fi
done

# Multiqc
$multiqc -f $workspace -o $workspace

# Upload the html file to s3
aws s3 cp $workspace/multiqc_report.html $output_address/$project_name/multiqc_report.html

# Download the html file to local--for display in the notebook
aws s3 cp $workspace/multiqc_report.html $notebook_dir
