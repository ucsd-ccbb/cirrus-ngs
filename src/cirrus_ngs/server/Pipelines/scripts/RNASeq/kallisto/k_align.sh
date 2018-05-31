#!/bin/bash
# This shell script uses kallisto alignment.

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
num_threads=${11}

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'k_align.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi
    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/ --quiet
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##


##ALIGNMENT##
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$kallisto/kallisto quant -i $kallisto_index -o $workspace/ -t $num_threads --single -l 50 -s 20 \
        $workspace/$fastq_end1$file_suffix" $JOB_NAME $status_file | tee $workspace/$fastq_end1.txt
else
    check_exit_status "$kallisto/kallisto quant -i $kallisto_index -o $workspace/ -t $num_threads \
        $workspace/$fastq_end1$file_suffix $workspace/$fastq_end2$file_suffix" $JOB_NAME $status_file | tee $workspace/$fastq_end1.txt
fi
##END_ALIGNMENT##

## Rename files ##
mv $workspace/abundance.h5 $workspace/$fastq_end1.abundance.h5
mv $workspace/abundance.tsv $workspace/$fastq_end1.abundance.tsv
mv $workspace/run_info.json $workspace/$fastq_end1.run_info.json

##UPLOAD##
aws s3 cp $workspace $output_address/ --exclude "*"  --include "$fastq_end1.abundance.*" --include "$fastq_end1.run_info.json" --include "$fastq_end1.txt" --recursive --quiet
##END_UPLOAD##
