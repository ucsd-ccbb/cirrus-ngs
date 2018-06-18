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
chromosome=${12}

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/"haplo.$chromosome.log"
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


##HAPLOTYPE##
check_exit_status "$bedtools genomecov -split -ibam $workspace/$fastq_end1$file_suffix \
    -bga -g $genome_fai -max 70001 > $workspace/$fastq_end1.final.$chromosome.bed" $JOB_NAME"_$chromosome" $status_file

check_exit_status "$java -Xms454m -Xmx8g -XX:+UseSerialGC -Djava.io.tmpdir=$workspace/temp \
    -jar $gatk -T HaplotypeCaller -R $genome_fasta \
    -I $workspace/$fastq_end1$file_suffix \
    -L $workspace/$fastq_end1.final.$chromosome.bed \
    --out $workspace/$fastq_end1.$chromosome.g.vcf \
    -nct $num_threads \
    -G none \
    --annotation BaseQualityRankSumTest \
    --annotation FisherStrand \
    --annotation GCContent \
    --annotation HomopolymerRun \
    --annotation MappingQualityRankSumTest \
    --annotation MappingQualityZero \
    --annotation QualByDepth \
    --annotation ReadPosRankSumTest \
    --annotation RMSMappingQuality \
    --annotation DepthPerAlleleBySample \
    --annotation Coverage \
    --annotation ClippingRankSumTest \
    --standard_min_confidence_threshold_for_calling 30.0 \
    -ERC GVCF \
    --dbsnp $dbsnp" $JOB_NAME"_$chromosome" $status_file

check_exit_status "check_outputs_exist $workspace/$fastq_end1.$chromosome.g.vcf" $JOB_NAME"_$chromosome" $status_file
##END_HAPLOTYPE##

##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.$chromosome.g.vcf" --recursive
##END_UPLOAD##

rm -r $workspace
