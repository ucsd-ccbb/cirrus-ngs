#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
chip_sample=$4
input_sample=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped

#logging
log_dir=$log_dir/$chip_sample
mkdir -p $log_dir
log_file=$log_dir/'findpeaks.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$chip_sample
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

pair_base_name=$chip_sample
rm -r $workspace/tags_$chip_sample $workspace/tags_$input_sample 

##DOWNLOAD##
if [ ! -d $workspace/tags_$chip_sample ] || [ "$input_sample" != "NULL" ]
then
    mkdir -p $workspace/tags_$chip_sample
    mkdir -p $workspace/tags_$input_sample

    #always download forward reads
    aws s3 cp $input_address/$project_name/$chip_sample/tags_$chip_sample $workspace/tags_$chip_sample --recursive

    #used to check if input_sample is the reverse reads or actually the input_sample
    aws s3 ls $input_address/$project_name/$input_sample/tags_$input_sample &>/dev/null

    if [ $? -eq 0 ]
    then
        pair_base_name=$chip_sample'_vs_'$input_sample
        aws s3 cp $input_address/$project_name/$input_sample/tags_$input_sample $workspace/tags_$input_sample --recursive
    fi
fi
##END_DOWNLOAD##


##FINDPEAKS## 
if [ "$pair_base_name" == "$chip_sample'_vs_'$input_sample" ]
then
    $find_peaks $workspace/tags_$chip_sample -style $style -o $workspace/$pair_base_name$style_ext -i $workspace/tags_$input_sample
else
    $find_peaks $workspace/tags_$chip_sample -style $style -o $workspace/$pair_base_name$style_ext
fi
##END_FINDPEAKS##


##UPLOAD##
aws s3 cp $workspace/$pair_base_name$style_ext $output_address/
##END_UPLOAD##
