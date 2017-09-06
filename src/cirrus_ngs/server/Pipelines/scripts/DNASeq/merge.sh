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
num_threads=${11}
chromosome_list=${12}
do_vcf_merging=${13} #either "True" or "False"

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'merge.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace
mkdir -p $workspace/temp

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

if [ "$do_vcf_merging" == "False" ]
then
    JOB_NAME=$JOB_NAME"_bamonly"
fi
check_step_already_done $JOB_NAME $status_file 

##DOWNLOAD##
downloads_needed="False"
for chrom in $chromosome_list
do
    if [ ! -f $workspace/$fastq_end1.final.$chrom.bam ] || [ ! -f $workspace/$fastq_end1.$chrom$file_suffix ]
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
        aws s3 cp $input_address/$fastq_end1.final.$chrom.bam $workspace/
        if [ "$do_vcf_merging" == "True" ]
        then
            aws s3 cp $input_address/$fastq_end1.$chrom$file_suffix $workspace/
        fi
    done
fi
##END_DOWNLOAD##

bam_file_list=""
vcf_file_list=""

for chrom in $chromosome_list
do
    bam_file_list=$bam_file_list"$workspace/$fastq_end1.final.$chrom.bam "
    vcf_file_list=$vcf_file_list"$workspace/$fastq_end1.$chrom$file_suffix "
done


##MERGE##
check_exit_status "$sambamba merge -t $num_threads $workspace/$fastq_end1.final.bam $bam_file_list" $JOB_NAME $status_file

if [ "$do_vcf_merging" == "True" ]
then
    check_exit_status "$vcf_concat $vcf_file_list > $workspace/$fastq_end1.raw.vcf" $JOB_NAME $status_file
    check_exit_status "$vcf_sort $workspace/temp $workspace/$fastq_end1.raw.vcf > $workspace/$fastq_end1.merged.vcf" $JOB_NAME $status_file
fi
##END_MERGE##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.bam" \
    --include "$fastq_end1.merged.vcf" \
    --recursive
##END_UPLOAD##
