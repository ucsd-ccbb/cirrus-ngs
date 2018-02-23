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
num_threads=${11}    #number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'gatk_align.log'
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

    #always download forward reads
    check_exit_status "aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/" $JOB_NAME $status_file

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        check_exit_status "aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/" $JOB_NAME $status_file
    fi
fi
##END_DOWNLOAD##

#TODO check human_genome, i'm pretty sure it's empty
# Star align
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
        --twopassMode Basic --readFilesIn $workspace/$fastq_end1$file_suffix \
        --outFileNamePrefix $workspace/$fastq_end1." $JOB_NAME $status_file

else
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
        --twopassMode Basic --readFilesIn $workspace/$fastq_end1$file_suffix $workspace/$fastq_end2$file_suffix \
        --outFileNamePrefix $workspace/$fastq_end1." $JOB_NAME $status_file
fi

check_exit_status "$samtools view -Sb $workspace/$fastq_end1.Aligned.out.sam > \
    $workspace/$fastq_end1.Aligned.out.bam" $JOB_NAME $status_file

#TODO: ask about sambamba
check_exit_status "$samtools sort -m 2G -@ 4 $workspace/$fastq_end1.Aligned.out.bam \
    $workspace/$fastq_end1.Aligned.out.sorted" $JOB_NAME $status_file

check_exit_status "$samtools index $workspace/$fastq_end1.Aligned.out.sorted.bam" $JOB_NAME $status_file

check_exit_status "$samtools stats $workspace/$fastq_end1.Aligned.out.sorted.bam > $workspace/$fastq_end1.txt" $JOB_NAME $status_file
# End star align



# Upload
aws s3 cp $workspace $output_address/ --exclude "*" --include "$fastq_end1.Aligned.out.sorted.bam*" --include "$fastq_end1.txt" --recursive
##END_UPLOAD##
