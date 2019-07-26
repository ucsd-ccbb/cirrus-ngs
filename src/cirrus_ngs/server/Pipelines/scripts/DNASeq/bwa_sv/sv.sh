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

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/"sv.log"
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

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/ --quiet
fi

if [ ! -f $workspace/$fastq_end1$file_suffix.bai ]
then
    check_exit_status "$sambamba index -t $num_threads $workspace/$fastq_end1$file_suffix \
      $workspace/$fastq_end1$file_suffix.bai" $JOB_NAME $status_file
fi

##END_DOWNLOAD##
##MANTA##
echo "running Manta..."
$python2 $manta --bam $workspace/$fastq_end1$file_suffix --referenceFasta $genome_fasta --runDir $workspace/Manta
$python2 $workspace/Manta/runWorkflow.py       -m local -j 16

##Run Lumpy##
# Extract the discordant paired-end alignments.
$samtools view -b -F 1294 $workspace/$fastq_end1$file_suffix > $workspace/$fastq_end1".discordants.unsorted"$file_suffix

# Extract the split-read alignments
$samtools view -h $workspace/$fastq_end1$file_suffix | $python2 $lumpy_extract -i stdin | $samtools view -Sb - > $workspace/$fastq_end1".splitters.unsorted"$file_suffix

# Sort both alignments
$samtools sort $workspace/$fastq_end1".discordants.unsorted"$file_suffix $workspace/$fastq_end1".discordants" &
$samtools sort $workspace/$fastq_end1".splitters.unsorted"$file_suffix $workspace/$fastq_end1".splitters" &
wait

# exclude low complexity regions
$lumpyexpress \
    -B $workspace/$fastq_end1$file_suffix \
    -S $workspace/$fastq_end1".splitters.bam" \
    -D $workspace/$fastq_end1".discordants.bam" \
    -x $exclude_bed \
    -o $workspace/$fastq_end1".lumpy.vcf"

# $bgzip $workspace/$fastq_end1".lumpy.vcf"
$sort_vcf $workspace/$fastq_end1".lumpy.vcf" | $bgzip > $workspace/$fastq_end1".lumpy.vcf.gz"
$tabix -p vcf $workspace/$fastq_end1".lumpy.vcf.gz"

# Concatenate SNV vcf
$concat_vcf $workspace/$fastq_end1".*.vcf.gz" > $workspace/$fastq_end1".snv.vcf"
$sort_vcf $workspace/$fastq_end1".snv.vcf" | $bgzip > $workspace/$fastq_end1".snv.vcf.gz"
$tabix -p vcf $workspace/$fastq_end1".snv.vcf.gz"

# Create PED file
echo -e "0\t$fastq_end1\t0\t0\t1" > $workspace/$fastq_end1".ped"

# Run SV2

$sv2 \
        -g mm10 \
        -i $workspace/$fastq_end1$file_suffix \
        -v $workspace/Manta/results/variants/diploidSV.vcf.gz $workspace/$fastq_end1".lumpy.vcf.gz" \
        -snv $workspace/$fastq_end1".vcf.gz" \
        -ped $workspace/$fastq_end1".ped" \
        -O $workspace

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $workspace/sv2_genotypes/sv2_genotypes.vcf > $workspace/sv2_genotypes/sv2_genotypes.filt.vcf

##UPLOAD##
aws s3 cp $workspace/Manta $output_address/Manta/ --recursive --quiet
aws s3 cp $workspace/sv2_features $output_address/sv2/sv2_features --recursive --quiet
aws s3 cp $workspace/sv2_genotypes $output_address/sv2/sv2_genotypes --recursive --quiet
aws s3 cp $workspace/sv2_preprocessing $output_address/sv2/sv2_preprocessing --recursive --quiet
##END_UPLOAD##

rm -r $workspace

