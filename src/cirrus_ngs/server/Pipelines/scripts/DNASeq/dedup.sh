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

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'dedup.log'
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
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet
    gunzip -q $workspace/$fastq_end1$download_suffix
fi
##END_DOWNLOAD##


##MARKDUPLICATES##
check_exit_status "java -jar -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g $picard_mark_duplicates \
    INPUT=$workspace/$fastq_end1$file_suffix OUTPUT=$workspace/$fastq_end1.dedup.bam \
    METRICS_FILE=$workspace/$fastq_end1.matrics.txt AS=true \
    VALIDATION_STRINGENCY=LENIENT" $JOB_NAME $status_file

check_exit_status "$sambamba index -t $num_threads $workspace/$fastq_end1.dedup.bam \
    $workspace/$fastq_end1.dedup.bam.bai" $JOB_NAME $status_file

check_exit_status "check_outputs_exist $workspace/$fastq_end1.dedup.bam \
    $workspace/$fastq_end1.dedup.bam.bai $workspace/$fastq_end1.matrics.txt" $JOB_NAME $status_file
##END_MARKDUPLICATES##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.dedup.bam*" --include "$fastq_end1.matrics.txt" --recursive --quiet
##END_UPLOAD##
