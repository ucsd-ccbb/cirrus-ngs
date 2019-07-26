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
chromosome=${13}

#logging
log_dir=$log_dir/$group_name
mkdir -p $log_dir
log_file=$log_dir/"filter.$chromosome.log"
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$group_name"_"$chromosome
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME"_$chromosome" $status_file

##DOWNLOAD##
if [ ! -f $workspace/$group_name$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and bam files
    aws s3 cp $input_address/$group_name/$group_name.$chromosome$download_suffix $workspace/ --quiet
    aws s3 cp $input_address/$group_name/$group_name.$chromosome$download_suffix".tbi" $workspace/ --quiet
fi
##END_DOWNLOAD##


##VARIANTFILTERING##
if [ -z $hapmap ] || [ -z $G1000_snps ] || [ -z $omni ]
then
    check_exit_status "$java -Xmx8g -jar $gatk \
        VariantRecalibrator \
        -R $genome_fasta \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_indels \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an FS -an SOR -an DP \
        -mode BOTH \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -V $workspace/$group_name.$chromosome.vcf.gz \
        -O $workspace/$group_name.snp.indel.$chromosome.vcf.gz \
        --tranches-file $workspace/$group_name.snp.indel.$chromosome.tranches" $JOB_NAME"_$chromosome" $status_file
else
    check_exit_status "$java -Xmx8g -jar $gatk \
        VariantRecalibrator \
        -R $genome_fasta \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_snps \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_indels \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
        -an QD -an MQ -an FS -an SOR -an DP \
        -mode BOTH \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -V $workspace/$group_name.$chromosome.vcf.gz \
        -O $workspace/$group_name.snp.indel.$chromosome.vcf.gz \
        --tranches-file $workspace/$group_name.snp.indel.$chromosome.tranches" $JOB_NAME"_$chromosome" $status_file
fi

check_exit_status "$java -Xmx8g -jar $gatk \
    ApplyVQSR \
    -V $workspace/$group_name.$chromosome.vcf.gz \
    -O $workspace/$group_name.snp.indel.vqsr.$chromosome.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $workspace/$group_name.snp.indel.$chromosome.tranches \
    --recal-file $workspace/$group_name.snp.indel.$chromosome.vcf.gz \
    -mode BOTH" $JOB_NAME"_$chromosome" $status_file


##END_VARIANTFILTERING##

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "`ls $workspace`"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

rm $workspace/$group_name.$chromosome.vcf.gz
rm $workspace/$group_name.$chromosome.vcf.gz.tbi

##UPLOAD##
aws s3 cp $workspace $output_address --recursive --quiet
##END_UPLOAD##

##CLEAN##
rm -r $workspace
##END_CLEAN##
