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

#logging
mkdir -p $log_dir
log_file=$log_dir/'k_diff_cal.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow
mkdir -p $workspace

echo $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

# Download group text file
if [ ! -f $workspace/group.txt ] || [ ! -f $workspace/all_gene_counts.txt ] 
then

    aws s3 cp $input_address/group.txt $workspace/group.txt
    aws s3 cp $input_address/all_gene_counts.txt $workspace/all_gene_counts.txt
fi
##END_DOWNLOAD##

# Call R script
check_exit_status "Rscript /shared/workspace/Pipelines/scripts/RNASeq/kallisto/RNA-seq_limma.R \
    $workspace/group.txt \
    $workspace/all_gene_counts.txt $workspace/" $JOB_NAME $status_file

echo "##########################"
ls $workspace/
echo "##########################"
