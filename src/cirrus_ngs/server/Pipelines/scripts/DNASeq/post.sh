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
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME"_$chromosome" $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ] || [ ! -f $workspace/$fastq_end1$file_suffix.bai ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/
    aws s3 cp $input_address/$fastq_end1$download_suffix.bai $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix
fi
##END_DOWNLOAD##


##POSTALIGN##
#picard - mark duplicates
check_exit_status "$java -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g -jar $mark_duplicates \
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
if [ -z $mills ] || [ -z $indels_1000G ] || [ -z $hapmap ]
then
    rtc_known="--known $indels --known $dbsnp"
    ir_known="-known $indels -known $dbsnp"
    br_known="--knownSites $indels --knownSites $dbsnp"
else
    rtc_known="--known $indels_1000G --known $mills"
    ir_known="-known $indels_1000G -known $mills"
    br_known="--knownSites $dbsnp --knownSites $hapmap"
fi

#gatk - realigner target creator
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx4g \
    $gatk -T RealignerTargetCreator -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    $rtc_known \
    -o $workspace/$fastq_end1.$chromosome.interval_list -L $chromosome \
    --allow_potentially_misencoded_quality_scores" $JOB_NAME"_$chromosome" $status_file

#gatk - indel realigner
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx4g \
    $gatk -T IndelRealigner -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    $ir_known \
    --targetIntervals $workspace/$fastq_end1.$chromosome.interval_list \
    -o $workspace/$fastq_end1.realign.$chromosome.bam -L $chromosome \
    --allow_potentially_misencoded_quality_scores" $JOB_NAME"_$chromosome" $status_file

#gatk - base recalibrator
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk -T BaseRecalibrator -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    $br_known \
    -o $workspace/$fastq_end1.$chromosome.grp \
    -dcov 2 -L $workspace/$fastq_end1.$chromosome.dedup.bed \
    --allow_potentially_misencoded_quality_scores" $JOB_NAME"_$chromosome" $status_file

#gatk - print reads
check_exit_status "$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk -T PrintReads -R $genome_fasta \
    -I $workspace/$fastq_end1.realign.$chromosome.bam \
    --BQSR $workspace/$fastq_end1.$chromosome.grp \
    -o $workspace/$fastq_end1.final.$chromosome.bam \
    -rf BadCigar -L $chromosome \
    --allow_potentially_misencoded_quality_scores" $JOB_NAME"_$chromosome" $status_file

check_exit_status "$sambamba index  $workspace/$fastq_end1.final.$chromosome.bam \
    $workspace/$fastq_end1.final.$chromosome.bam.bai" $JOB_NAME"_$chromosome" $status_file

##END_POSTALIGN##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.$chromosome.bam*" --recursive
##END_UPLOAD##
