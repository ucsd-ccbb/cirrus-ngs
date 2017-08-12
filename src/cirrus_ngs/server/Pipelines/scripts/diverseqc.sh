#!/bin/bash

project_name=$1
root_dir=$2            # /scratch
input_address=$3    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$4   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$5
array=( $@ )
len=${#array[@]}
multiqc_file_list=${array[@]:5:$len}  # multiqc files: unknown number of args

mkdir -p $log_dir
log_file=$log_dir/'diverseqc.log'
exec 1>>$log_file
exec 2>>$log_file

workspace=$root_dir/$project_name
software=/shared/workspace/software
multiqc=$software/MultiQC/multiqc

mkdir -p $workspace

# Download files to perform multiqc on from s3
for i in ${multiqc_file_list//,/ }; do
    echo "In shell: "$i
    if [ ! -f $workspace/$i ]
    then
        aws s3 cp $input_address/$i $workspace
    fi
done

# Multiqc
# -f will overwrite the old report with the same name
$multiqc -f $workspace -o $workspace

# Upload the html file to s3
aws s3 cp $workspace/multiqc_report.html $output_address/$project_name/multiqc_report.html


