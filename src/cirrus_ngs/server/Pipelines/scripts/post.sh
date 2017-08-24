#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
numb_threads=${10}
chromosome=${11}

#logging
mkdir -p $log_dir
log_file=$log_dir/'post.log'
exec 1>>$log_file
exec 2>>$log_file

. /shared/workspace/software/software.conf

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

#reference files
genome_fai=$software_dir/sequences/Hsapiens/ucsc.hg19.fasta.fai
genome_fasta=$software_dir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp=$software_dir/variation/dbsnp_138.hg19.vcf
mills=$software_dir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
G1000=$software_dir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf
hapmap=$software_dir/variation/hapmap_3.3.hg19.sites.vcf

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
$java -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g -jar $mark_duplicates \
    INPUT=$workspace/$fastq_end1$file_suffix \
    OUTPUT=$workspace/$fastq_end1.$chromosome.dedup.bam \
    METRICS_FILE=$workspace/$fastq_end1.$chromosome.metrics.txt \
    AS=true VALIDATION_STRINGENCY=LENIENT

#sambamba - indexing
$sambamba index  $workspace/$fastq_end1.$chromosome.dedup.bam \
    $workspace/$fastq_end1.$chromosome.dedup.bam.bai

#bedtools - create bed file
$bedtools genomecov -split -ibam $workspace/$fastq_end1.$chromosome.dedup.bam \
    -bga -g $genome_fai -max 70001 > $workspace/$fastq_end1.$chromosome.dedup.bed

#gatk - realigner target creator
$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx4g \
    $gatk -T RealignerTargetCreator -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    --known $G1000 --known $mills \
    -o $workspace/$fastq_end1.$chromosome.interval_list -L chr$chromosome \
    --allow_potentially_misencoded_quality_scores

#gatk - indel realigner
$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx4g \
    $gatk -T IndelRealigner -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    -known $G1000 -known $mills \
    --targetIntervals $workspace/$fastq_end1.$chromosome.interval_list \
    -o $workspace/$fastq_end1.realign.$chromosome.bam -L chr$chromosome \
    --allow_potentially_misencoded_quality_scores

#gatk - base recalibrator
$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk -T BaseRecalibrator -R $genome_fasta \
    -I $workspace/$fastq_end1.$chromosome.dedup.bam \
    --knownSites $dbsnp --knownSites $hapmap \
    -o $workspace/$fastq_end1.$chromosome.grp \
    -dcov 2 -L $workspace/$fastq_end1.$chromosome.dedup.bed \
    --allow_potentially_misencoded_quality_scores

#gatk - print reads
$java -d64 -jar -Djava.io.tmpdir=$workspace/temp -Xmx8g \
    $gatk -T PrintReads -R $genome_fasta \
    -I $workspace/$fastq_end1.realign.$chromosome.bam \
    --BQSR $workspace/$fastq_end1.$chromosome.grp \
    -o $workspace/$fastq_end1.final.$chromosome.bam \
    -rf BadCigar -L chr$chromosome \
    --allow_potentially_misencoded_quality_scores
    #--no_pg_tag

$sambamba index  $workspace/$fastq_end1.final.$chromosome.bam \
    $workspace/$fastq_end1.final.$chromosome.bam.bai

##END_POSTALIGN##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.$chromosome.bam*" --recursive
##END_UPLOAD##
