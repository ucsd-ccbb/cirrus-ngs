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
num_threads=${10}
chromosome_list=${11}

#logging
mkdir -p $log_dir
log_file=$log_dir/'merge.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
software_dir=/shared/workspace/software
sambamba=$software_dir/sambamba/0.4.7/bin/sambamba
java=$software_dir/java/jre1.8.0_144/bin/java
gatk=$software_dir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar

genome_fasta=$software_dir/sequences/Hsapiens/ucsc.hg19.fasta

vcf_concat=$software_dir/vcftools_0.1.12b/bin/vcf-concat
vcf_sort=$software_dir/vcftools_0.1.12b/bin/vcf-sort
bgzip=$software_dir/tabix-0.2.6/bgzip
tabix=$software_dir/tabix-0.2.6/tabix
mkdir -p $workspace


##DOWNLOAD##
downloads_needed="False"
for chrom in $chromosome_list
do
    if [ ! -f $workspace/$fastq_end1.final.$chrom.bam ] || [ ! -f $workspace/$fastq_end1.$chrom.g.vcf.gz ] || [ ! -f $workspace/$fastq_end1.$chrom.g.vcf.gz.tbi ]
    then
        downloads_needed="True"
    fi
done

if [ "$downloads_needed" == "True" ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and bam files
    for chrom in $chromosome_list
    do
        aws s3 cp $input_address/$fastq_end1.final.$chrom.bam $workspace/
        aws s3 cp $input_address/$fastq_end1.$chrom.g.vcf.gz $workspace/
        aws s3 cp $input_address/$fastq_end1.$chrom.g.vcf.gz.tbi $workspace/
    done
fi
##END_DOWNLOAD##

bam_file_list=""
variant_list=""

for chrom in $chromosome_list
do
    bam_file_list=$bam_file_list"$workspace/$fastq_end1.final.$chrom.bam "
    variant_list=$variant_list"--variant $workspace/$fastq_end1.$chrom.g.vcf.gz "
done


##MERGE##
$sambamba merge -t $num_threads $workspace/$fastq_end1.final.bam $bam_file_list 

$java -Xmx2g -Djava.io.tmpdir=$workspace/temp \
    -jar $gatk \
    -T CombineGVCFs \
    -R $genome_fasta \
    $variant_list \
    -o $workspace/$fastq_end1.merged.vcf.gz

##END_MERGE##

##UPLOAD##
#--include "$fastq_end1.raw.vcf.gz" \
aws s3 cp $workspace $output_address --exclude "*" --include "$fastq_end1.final.bam" \
    --include "$fastq_end1.merged.vcf*" \
    --recursive
##END_UPLOAD##
