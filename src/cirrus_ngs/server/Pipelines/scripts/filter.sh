#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
group_name=$4       ##NEW## name of group being combined
fastq_end2=$5       ##NEW## always null now
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
files_in_group=${10}    #all files in current group
num_threads=${11}

output_address=$output_address/$project_name/$group_name

#logging
mkdir -p $log_dir
log_file=$log_dir/'filter_variant.log'
exec 1>>$log_file
exec 2>>$log_file

#prepare output directories
workspace=$root_dir/$project_name/$group_name
software_dir=/shared/workspace/software
java=$software_dir/java/jre1.8.0_144/bin/java
gatk=$software_dir/gatk/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
bgzip=$software_dir/tabix-0.2.6/bgzip
tabix=$software_dir/tabix-0.2.6/tabix
mkdir -p $workspace

#reference files
genome_fai=$software_dir/sequences/Hsapiens/ucsc.hg19.fasta.fai
genome_fasta=$software_dir/sequences/Hsapiens/ucsc.hg19.fasta
dbsnp=$software_dir/variation/dbsnp_138.hg19.vcf
mills=$software_dir/variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
omni=$software_dir/variation/1000G_omni2.5.hg19.sites.vcf
hapmap=$software_dir/variation/hapmap_3.3.hg19.sites.vcf
G1000snps=$software_dir/variation/1000G_phase1.snps.high_confidence.hg19.sites.vcf
G1000indels=$software_dir/variation/1000G_phase1.indels.hg19.sites.vcf


##DOWNLOAD##
if [ ! -f $workspace/$group_name$file_suffix ] || [ ! -f $workspace/$group_name$file_suffix.idx ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and bam files
    aws s3 cp $input_address/$group_name$file_suffix $workspace/
    aws s3 cp $input_address/$group_name$file_suffix.tbi $workspace/
fi
##END_DOWNLOAD##


##VARIANTFILTERING##
$java -Xmx4g -jar $gatk \
    -nt $num_threads \
    -R $genome_fasta -l debug \
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
    --disable_auto_index_creation_and_locking_when_reading_rods

$java -Xmx4g -jar $gatk \
    -nt $num_threads \
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
    --disable_auto_index_creation_and_locking_when_reading_rods

$java -Xmx4g -jar $gatk \
    -nt $num_threads \
    -R $genome_fasta \
    -T ApplyRecalibration \
    -mode SNP \
    --ts_filter_level 99.0 \
    -input $workspace/$group_name$file_suffix \
    -recalFile $workspace/$group_name.snp.recal \
    -tranchesFile $workspace/$group_name.snp.tranches \
    -o $workspace/$group_name.snpAr.vcf

$java -Xmx4g -jar $gatk \
    -nt $num_threads \
    -R $genome_fasta \
    -T ApplyRecalibration \
    -mode indel \
    --ts_filter_level 99.0 \
    -input $workspace/$group_name.snpAr.vcf \
    -recalFile $workspace/$group_name.indel.recal \
    -tranchesFile $workspace/$group_name.indel.tranches \
    -o $workspace/$group_name.snpAr.indelAr.vcf

$java -Xmx4g -jar $gatk \
    -nt $num_threads \
    -R $genome_fasta \
    -T SelectVariants \
    --excludeNonVariants \
    --excludeFiltered \
    --variant $workspace/$group_name.snpAr.indelAr.vcf \
    --out $workspace/$group_name.vqsr.vcf

##END_VARIANTFILTERING##

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "`ls $workspace`"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
