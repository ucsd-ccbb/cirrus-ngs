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
all_samples=${11}

#logging
log_dir=$log_dir
mkdir -p $log_dir
log_file=$log_dir/'anno_filt.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow
tempDir=$workspace/tmp

mkdir -p $tempDir

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

##DOWNLOAD##
check_exit_status "aws s3 cp $input_address/$project_name'_merged_alu'$file_suffix'.gz' $workspace/ --quiet" $JOB_NAME $status_file
check_exit_status "aws s3 cp $input_address/$project_name'_merged_nonalu'$file_suffix'.gz' $workspace/ --quiet" $JOB_NAME $status_file
##END_DOWNLOAD##

check_exit_status "$oncotator -i VCF --skip-no-alt --db-dir $oncotator_db $workspace/$project_name'_merged_alu'$file_suffix'.gz' $workspace/$project_name'_merged_alu'.maf hg19" $JOB_NAME $status_file

check_exit_status "$oncotator -i VCF --skip-no-alt --db-dir $oncotator_db $workspace/$project_name'_merged_nonalu'$file_suffix'.gz' $workspace/$project_name'_merged_nonalu'.maf hg19" $JOB_NAME $status_file

$Rscript $annotation_filt $workspace/$project_name'_merged_anno_filt.maf' $workspace/$project_name'_merged_alu.maf' $workspace/$project_name'_merged_nonalu.maf'

check_exit_status "check_outputs_exist $workspace/$project_name'_merged_anno_filt.maf'" $JOB_NAME $status_file

##UPLOAD##
aws s3 cp $workspace/$project_name'_merged_anno_filt.maf' $output_address/ --quiet
aws s3 cp $workspace/$project_name'_merged_alu.maf' $output_address/ --quiet
aws s3 cp $workspace/$project_name'_merged_nonalu.maf' $output_address/ --quiet
##END_UPLOAD##

##CLEAN##
rm -r $workspace
##END_CLEAN##
