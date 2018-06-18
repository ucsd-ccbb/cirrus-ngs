#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
chip_sample=$5
input_sample=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped

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

rm -r $workspace/tags_$chip_sample $workspace/tags_$input_sample &>/dev/null

pair_base_name=$chip_sample

##DOWNLOAD##
if [ ! -d $workspace/tags_$chip_sample ] || [ "$pairs_exist" == "True" ]
then
    mkdir -p $workspace/tags_$chip_sample

    if [ "$pairs_exist" == "True" ]
    then
        mkdir -p $workspace/tags_$input_sample
        pair_base_name=$chip_sample'_vs_'$input_sample
        aws s3 cp $input_address/$input_sample/tags_$input_sample $workspace/tags_$input_sample --recursive
        input_address=$input_address/$chip_sample
    fi

    aws s3 cp $input_address/tags_$chip_sample $workspace/tags_$chip_sample --recursive
fi
##END_DOWNLOAD##


##FINDPEAKS##
if [ "$pairs_exist" == "True" ]
then
    check_exit_status "$find_peaks $workspace/tags_$chip_sample -style $style \
        -o $workspace/$pair_base_name$style_ext -i $workspace/tags_$input_sample" $JOB_NAME $status_file
else
    check_exit_status "$find_peaks $workspace/tags_$chip_sample -style $style \
        -o $workspace/$pair_base_name$style_ext" $JOB_NAME $status_file
fi

check_exit_status "check_outputs_exist $workspace/$pair_base_name$style_ext" $JOB_NAME $status_file
##END_FINDPEAKS##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$pair_base_name$style_ext" --recursive
##END_UPLOAD##

rm -r $workspace
