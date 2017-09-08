#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
group_name=$5
fastq_end2=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
files_in_group=${11}    #all files in current group
num_threads=${12}

#logging
log_dir=$log_dir/$group_name
mkdir -p $log_dir
log_file=$log_dir/'filter_variant.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$group_name
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$group_name$file_suffix ] || [ ! -f $workspace/$group_name$file_suffix.tbi ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and bam files
    aws s3 cp $input_address/$group_name/$group_name$file_suffix $workspace/
    aws s3 cp $input_address/$group_name/$group_name$file_suffix.tbi $workspace/
fi
##END_DOWNLOAD##


##VARIANTFILTERING##
check_exit_status "$java -Xmx8g -jar $gatk \
    -nt $num_threads \
    -R $genome_fasta \
    -T VariantRecalibrator \
    --maxGaussians 4 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000snps \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -input $workspace/$group_name.g.vcf.gz \
    -recalFile $workspace/$group_name.snp.recal \
    -tranchesFile $workspace/$group_name.snp.tranches \
    -rscriptFile $workspace/$group_name.snp.plots.R \
    --disable_auto_index_creation_and_locking_when_reading_rods" $JOB_NAME $status_file

check_exit_status "$java -Xmx8g -jar $gatk \
    -R $genome_fasta \
    -T VariantRecalibrator \
    --maxGaussians 4 \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000indels \
    -an DP -an FS -an ReadPosRankSum -an MQRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -input $workspace/$group_name$file_suffix \
    -recalFile $workspace/$group_name.indel.recal \
    -tranchesFile $workspace/$group_name.indel.tranches \
    -rscriptFile $workspace/$group_name.indel.plots.R \
    --disable_auto_index_creation_and_locking_when_reading_rods" $JOB_NAME $status_file

##END_VARIANTFILTERING##

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "`ls $workspace`"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
