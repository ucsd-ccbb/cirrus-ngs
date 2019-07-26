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
numb_threads=${11}
chromosome=${12}

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/"post.$chromosome.log"
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace/temp

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME"_$chromosome" $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ] || [ ! -f $workspace/$fastq_end1$file_suffix.bai ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet
    aws s3 cp $input_address/$fastq_end1$download_suffix.bai $workspace/ --quiet
fi
##END_DOWNLOAD##


##POSTALIGN##
#picard - mark duplicates
check_exit_status "$java -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g -jar $picard_mark_duplicates \
    INPUT=$workspace/$fastq_end1$file_suffix \
    OUTPUT=$workspace/$fastq_end1.$chromosome.dedup.bam \
    METRICS_FILE=$workspace/$fastq_end1.$chromosome.metrics.txt \
    AS=true VALIDATION_STRINGENCY=LENIENT" $JOB_NAME"_$chromosome" $status_file

#sambamba - indexing
check_exit_status "$sambamba index  $workspace/$fastq_end1.$chromosome.dedup.bam \
    $workspace/$fastq_end1.$chromosome.dedup.bam.bai" $JOB_NAME"_$chromosome" $status_file

#bedtools - create bed file
check_exit_status "$bedtools genomecov -split -ibam $workspace/$fastq_end1.$chromosome.dedup.bam \
    -bga -g $genome_fai -max 70001 > $workspace/$fastq_end1.$chromosome.dedup.bed" $JOB_NAME"_$chromosome" $status_file

#sets proper reference files depending on availablity
if [ -z $mills ] || [ -z $G1000_indels ] || [ -z $hapmap ]
then
    rtc_known="--known $dbsnp"
    ir_known="-known $dbsnp"
    br_known="--known-sites $dbsnp"
else
    rtc_known="--known $G1000_indels --known $mills"
    ir_known="-known $G1000_indels -known $mills"
    br_known="--known-sites $dbsnp --known-sites $hapmap"
fi

#gatk - base recalibrator
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk BaseRecalibrator -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    $br_known \
    -O $workspace/$fastq_end1.$chromosome.recalibration.table" $JOB_NAME"_$chromosome" $status_file

#gatk - apply bqsr
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk ApplyBQSR -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    --bqsr-recal-file $workspace/$fastq_end1.$chromosome.recalibration.table \
    -O $workspace/$fastq_end1.final.$chromosome.bam" $JOB_NAME"_$chromosome" $status_file

check_exit_status "$sambamba index  $workspace/$fastq_end1.final.$chromosome.bam \
    $workspace/$fastq_end1.final.$chromosome.bam.bai" $JOB_NAME"_$chromosome" $status_file

check_exit_status "check_outputs_exist $workspace/$fastq_end1.final.$chromosome.bam \
    $workspace/$fastq_end1.final.$chromosome.bam.bai" ${JOB_NAME}_$chromosome $status_file

##END_POSTALIGN##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.$chromosome.bam*" --recursive --quiet
##END_UPLOAD##

##CLEAN##
rm $workspace/$fastq_end1$file_suffix
rm $workspace/$fastq_end1$file_suffix.bai
rm $workspace/$fastq_end1.final.$chromosome.bam
rm $workspace/$fastq_end1.final.$chromosome.bam.bai
##END_CLEAN##
