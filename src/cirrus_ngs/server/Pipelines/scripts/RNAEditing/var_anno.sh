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
log_file=$log_dir/'var_anno.log'
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
for file in $all_samples; 
do
    	check_exit_status "aws s3 cp $input_address/$file/$file'_alu'$file_suffix $workspace/ --quiet" $JOB_NAME $status_file
    	check_exit_status "aws s3 cp $input_address/$file/$file'_nonalu'$file_suffix $workspace/ --quiet" $JOB_NAME $status_file
done
##END_DOWNLOAD##

for file in $all_samples;
do 
	$bgzip $workspace/$file'_alu'$file_suffix; $tabix -p vcf $workspace/$file'_alu'$file_suffix.gz;
   	$bgzip $workspace/$file'_nonalu'$file_suffix; $tabix -p vcf $workspace/$file'_nonalu'$file_suffix.gz; 
done

alu_vcf_list=""
nonalu_vcf_list=""

for file in $all_samples;
do
    	alu_vcf_list=$alu_vcf_list"$workspace/$file"_alu"$file_suffix.gz "
 	nonalu_vcf_list=$nonalu_vcf_list"$workspace/$file"_nonalu"$file_suffix.gz "
done

check_exit_status "$bcftools merge $alu_vcf_list -o $workspace/$project_name"_merged_alu"$file_suffix" $JOB_NAME $status_file
check_exit_status "$bcftools merge $nonalu_vcf_list -o $workspace/$project_name"_merged_nonalu"$file_suffix" $JOB_NAME $status_file

check_exit_status "$bgzip $workspace/$project_name"_merged_alu"$file_suffix; $tabix -p vcf $workspace/$project_name"_merged_alu"$file_suffix.gz;" $JOB_NAME $status_file
check_exit_status "$bgzip $workspace/$project_name"_merged_nonalu"$file_suffix; $tabix -p vcf $workspace/$project_name"_merged_nonalu"$file_suffix.gz;" $JOB_NAME $status_file

##UPLOAD## 
check_exit_status "aws s3 cp $workspace/$project_name'_merged_alu'$file_suffix'.gz' $output_address/ --quiet" $JOB_NAME $status_file
check_exit_status "aws s3 cp $workspace/$project_name'_merged_alu'$file_suffix'.gz.tbi' $output_address/ --quiet" $JOB_NAME $status_file
check_exit_status "aws s3 cp $workspace/$project_name'_merged_nonalu'$file_suffix'.gz' $output_address/ --quiet" $JOB_NAME $status_file
check_exit_status "aws s3 cp $workspace/$project_name'_merged_nonalu'$file_suffix'.gz.tbi' $output_address/ --quiet" $JOB_NAME $status_file
##END_UPLOAD##

rm -r $workspace

