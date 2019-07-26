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
log_file=$log_dir/'star_align.log'
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

#this is the suffix of the input from s3
download_suffix=$file_suffix

if [ "$is_zipped" == "True" ]
then
    download_suffix=$file_suffix".gz"
fi

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$download_suffix ]
then
    #always download forward reads
    check_exit_status "aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet" $JOB_NAME $status_file

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        check_exit_status "aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/ --quiet" $JOB_NAME $status_file
    fi
fi
##END_DOWNLOAD##

# Star align
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
        --twopassMode Basic --readFilesIn $workspace/$fastq_end1$download_suffix \
	--readFilesCommand zcat --quantMode GeneCounts --sjdbGTFfile $hg19_gtf \
	--alignIntronMax 200000 --alignMatesGapMax 200000 \
	--outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonical \
	--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $num_threads \
	--outSAMunmapped Within --outSAMattrRGline ID:"1234" SM:"$fastq_end1" PL:"ILLUMINA" LB:"TrueSeq" \
	--outSAMstrandField intronMotif --chimSegmentMin 25 --chimJunctionOverhangMin 25 \
	--limitBAMsortRAM 31311727203 --outFileNamePrefix $workspace/$fastq_end1." \
        $JOB_NAME $status_file

else
    check_exit_status "$STAR --runThreadN $num_threads --genomeDir $STAR_index \
        --twopassMode Basic --readFilesIn $workspace/$fastq_end1$download_suffix $workspace/$fastq_end2$download_suffix \
        --readFilesCommand zcat --quantMode GeneCounts --sjdbGTFfile $hg19_gtf \
        --alignIntronMax 200000 --alignMatesGapMax 200000 \
        --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $num_threads \
        --outSAMunmapped Within --outSAMattrRGline ID:"1234" SM:"$fastq_end1" PL:"ILLUMINA" LB:"TrueSeq" \
        --outSAMstrandField intronMotif --chimSegmentMin 25 --chimJunctionOverhangMin 25 \
        --limitBAMsortRAM 31311727203 --outFileNamePrefix $workspace/$fastq_end1." \
        $JOB_NAME $status_file

fi

check_exit_status "$samtools index $workspace/$fastq_end1.Aligned.sortedByCoord.out.bam" $JOB_NAME $status_file

check_exit_status "$samtools stats $workspace/$fastq_end1.Aligned.sortedByCoord.out.bam > $workspace/$fastq_end1.txt" $JOB_NAME $status_file
# End star align

# Upload
aws s3 cp $workspace $output_address/ --exclude "*" --include "$fastq_end1.Aligned.sortedByCoord.out.bam*" --include "$fastq_end1.txt" --recursive --quiet
##END_UPLOAD##

rm -r $workspace
