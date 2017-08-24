#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
num_threads=${10}
chromosome_list=${11}

#logging
mkdir -p $log_dir
log_file=$log_dir/'merge.log'
exec 1>>$log_file
exec 2>>$log_file

. /shared/workspace/software/software.conf

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace
mkdir -p $workspace/temp

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
        aws s3 cp $input_address/$fastq_end1.$chrom$file_suffix $workspace/
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
export PERL5LIB=/shared/workspace/software/vcftools_0.1.12b/perl/
$sambamba merge -t $num_threads $workspace/$fastq_end1.final.bam $bam_file_list 

$vcf_concat $vcf_file_list > $workspace/$fastq_end1.raw.vcf
$vcf_sort $workspace/temp $workspace/$fastq_end1.raw.vcf > $workspace/$fastq_end1.merged.vcf
##END_MERGE##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.bam" \
    --include "$fastq_end1.merged.vcf" \
    --recursive
##END_UPLOAD##
