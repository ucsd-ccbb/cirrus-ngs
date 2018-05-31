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
num_threads=${11}   # number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'cal_expression.log'
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
if [ ! -f $workspace/$fastq_end1$file_suffix ] || [ ! -f $workspace/$fastq_end2$file_suffix ] 
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/ --quiet
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##

# RSEM calculate expression
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$rsem --star --output-genome-bam --sort-bam-by-coordinate \
        --star-path $star_path -p $num_threads \
        $workspace/$fastq_end1$file_suffix \
        $rsem_index $workspace/$fastq_end1" $JOB_NAME $status_file

else
    # paired end
    check_exit_status "$rsem --star --output-genome-bam --sort-bam-by-coordinate \
        --star-path $star_path -p $num_threads \
        --paired-end $workspace/$fastq_end1$file_suffix \
        $workspace/$fastq_end2$file_suffix \
        $rsem_index $workspace/$fastq_end1" $JOB_NAME $status_file
fi

# perform samtools stats, for multiqc purposes
check_exit_status "$samtools stats $workspace/$fastq_end1.genome.bam > $workspace/$fastq_end1.txt" $JOB_NAME $status_file
check_exit_status "check_outputs_exist $workspace/$fastq_end1.txt $workspace/$fastq_end1.genes.results \
    $workspace/$fastq_end1.isoforms.results $workspace/$fastq_end1.stat \
    $workspace/$fastq_end1.transcript.sorted.bam $workspace/$fastq_end1.transcript.sorted.bam.bai \
    $workspace/$fastq_end1.genome.sorted.bam $workspace/$fastq_end1.genome.sorted.bam.bai" \
    $JOB_NAME $status_file


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.transcript.sorted.bam*" --include "$fastq_end1.txt" \
    --include "$fastq_end1.*.results" --include "$fastq_end1.genome.sorted.bam*" --recursive --quiet
aws s3 cp $workspace/$fastq_end1.stat $output_address/$fastq_end1.stat --recursive --quiet
##END_UPLOAD##
