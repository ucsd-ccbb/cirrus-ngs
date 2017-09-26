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
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'bowtie2.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace
basename=hairpin_human

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo $file_suffix
echo $is_zipped
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
    check_exit_status "aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/" $JOB_NAME $status_file
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        check_exit_status "aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/" $JOB_NAME $status_file
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi
fi
##END_DOWNLOAD##

##BOWTIE 2 ALIGNMENT##

# go to the directory of the index file
cd $bowtie2_index

# indexing a reference genome (already built, moved under software/bowtie2_index)
# cd ~
# $bowtie/bowtie2-build $hairpin_human_fa $basename

if [ "$fastq_end2" == "NULL" ]
then
    # single end
    check_exit_status "$bowtie2/bowtie2 --local -p 8 -q --phred33 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 \
     -x $basename -U $workspace/$fastq_end1$file_suffix \
    -S $workspace/$fastq_end1.sam" $JOB_NAME $status_file
else
    # paired end: using the same output file name ($fastq_end1.sam)
    check_exit_status "$bowtie2/bowtie2 --local -p 8 -q --phred33 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 \
    -x $basename -1 $workspace/$fastq_end1$file_suffix \
    -2 $workspace/$fastq_end2/$file_suffix -S $workspace/$fastq_end1.sam" $JOB_NAME $status_file
fi

##END BOWTIE 2 ##

# Produce text files to workspace, for multiqc analysis
check_exit_status "$samtools stats $workspace/$fastq_end1.sam > $workspace/$fastq_end1.txt" $JOB_NAME $status_file
echo "Finished samtools stats"

##UPLOAD##
# upload the sam files
aws s3 cp $workspace $output_address --exclude "*" --include "*.sam*" --recursive

# upload the txt files from samtool stats
aws s3 cp $workspace $output_address --exclude "*" --include "*.txt*" --recursive
##END_UPLOAD##

