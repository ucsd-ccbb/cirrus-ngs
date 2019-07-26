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
log_file=$log_dir/'align.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace/$fastq_end1

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

#this is the suffix of the input from s3
download_suffix=$file_suffix

if [ "$is_zipped" == "True" ]
then
    download_suffix=$file_suffix".gz"
fi

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ] || [ ! -f $workspace/$fastq_end2$file_suffix ] 
then
    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/ --quiet
    fi
fi
##END_DOWNLOAD##

# STAR alignment
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
	--outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  \
	--outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  \
	--alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  \
	--alignSJDBoverhangMin 1 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted \
	--quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted \
        --outFileNamePrefix $workspace/$fastq_end1/ --readFilesCommand zcat \
        --readFilesIn $workspace/$fastq_end1$download_suffix" $JOB_NAME $status_file
else
    # paired end
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
        --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD  \
        --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  \
        --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  \
        --alignSJDBoverhangMin 1 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM  --outSAMheaderHD \@HD VN:1.4 SO:unsorted \
	--outFileNamePrefix $workspace/$fastq_end1/ --readFilesCommand zcat \
        --readFilesIn $workspace/$fastq_end1$download_suffix $workspace/$fastq_end2$download_suffix" $JOB_NAME $status_file
fi

check_exit_status "check_outputs_exist $workspace/$fastq_end1/Aligned.toTranscriptome.out.bam \
	$workspace/$fastq_end1/Aligned.out.bam" $JOB_NAME $status_file

##UPLOAD##
aws s3 cp $workspace/$fastq_end1/ $output_address --recursive --quiet
##END_UPLOAD##

##CLEAN##
#rm -r $workspace
##END_CLEAN##
