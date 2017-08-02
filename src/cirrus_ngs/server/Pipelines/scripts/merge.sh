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

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software_dir=/shared/workspace/software
sambamba=$software_dir/sambamba/0.4.7/bin/sambamba
vcf_concat=$software_dir/vcftools_0.1.12b/bin/vcf-concat
vcf_sort=$software_dir/vcftools_0.1.12b/bin/vcf-sort
bgzip=$software_dir/tabix-0.2.6/bgzip
tabix=$software_dir/tabix-0.2.6/tabix
mkdir -p $workspace

##DOWNLOAD##
downloads_needed="False"
for chrom in $chromosome_list
do
    if [ ! -f $workspace/$fastq_end1.final.$chrom.bam ]
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
        aws s3 cp $input_address/$fastq_end1.final.$chrom.vcf.gz $workspace/
    done
fi
##END_DOWNLOAD##

bam_file_list=""
vcf_file_list=""

for chrom in $chromosome_list
do
    bam_file_list=$bam_file_list"$workspace/$fastq_end1.final.$chrom.bam "
    vcf_file_list=$vcf_file_list"$workspace/$fastq_end1.final.$chrom.vcf.gz "
done


##MERGE##
$sambamba merge -t $num_threads $workspace/$fastq_end1.final.bam $bam_file_list 
$vcf_concat $vcf_file_list | $bgzip -c > $workspace/$fastq_end1.final.vcf.gz
$tabix -p vcf $workspace/$fastq_end1.final.vcf.gz
##END_MERGE##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.bam" \
    --include "$fastq_end1.final.vcf.gz*" --include "$fastq_end1.raw.vcf.gz" \
    --recursive
##END_UPLOAD##
