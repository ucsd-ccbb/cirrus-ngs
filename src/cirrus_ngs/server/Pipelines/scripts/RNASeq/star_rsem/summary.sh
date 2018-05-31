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
#log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'summary.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

# Download files with extensions: .cnt, .genes.results, and .isoforms.results
for file in $all_samples
do
    aws s3 cp $input_address/$file/ $workspace/ --exclude "*" --include "*.cnt" --include "*.genes.results" --include "*.isoforms.results" --recursive --quiet
done

# Call the python parser files and upload the output file to s3
$python /shared/workspace/Pipelines/util/star_rsem/RSEM_count_parser.py $workspace
aws s3 cp $workspace $output_address --exclude "*" --include "*all_counts_results.txt*" --recursive --quiet

$python /shared/workspace/Pipelines/util/star_rsem/RSEM_gene_parser.py $workspace
aws s3 cp $workspace $output_address --exclude "*" --include "*all_genes_results.txt*" --recursive --quiet

$python /shared/workspace/Pipelines/util/star_rsem/RSEM_isoform_parser.py $workspace
aws s3 cp $workspace $output_address --exclude "*" --include "*all_isoforms_results.txt*" --recursive --quiet