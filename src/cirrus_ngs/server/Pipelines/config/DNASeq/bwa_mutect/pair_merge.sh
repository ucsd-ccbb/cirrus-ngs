#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
normal_sample=$5
tumor_sample=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
num_threads=${11}

#logging
log_dir=$log_dir/$normal_sample
mkdir -p $log_dir
log_file=$log_dir/'pair_merge.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$normal_sample
mkdir -p $workspace
mkdir -p $workspace/temp

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file 
pair_base_name=$normal_sample'_vs_'$tumor_sample

##DOWNLOAD##
downloads_needed="False"
for chrom in $chromosome_list
do
    if [ ! -f $workspace/$pair_base_name.$chrom$file_suffix ]
    then
        downloads_needed="True"
    fi
done

if [ "$downloads_needed" == "True" ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and bam files
    for chrom in $chromosome_list
    do
        aws s3 cp $input_address/$normal_sample/$pair_base_name.$chrom$file_suffix $workspace/
    done
fi
##END_DOWNLOAD##

vcf_file_list=""

for chrom in $chromosome_list
do
    vcf_file_list=$vcf_file_list"$workspace/$pair_base_name.$chrom$file_suffix "
done


##MERGE##
check_exit_status "$vcf_concat $vcf_file_list > $workspace/$pair_base_name.raw.vcf" $JOB_NAME $status_file
check_exit_status "$python $vcf_sort $workspace/$pair_base_name.raw.vcf '$chromosome_list' -o $workspace/$pair_base_name.merged.vcf" $JOB_NAME $status_file
##END_MERGE##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$pair_base_name.merged.vcf" --recursive
##END_UPLOAD##
