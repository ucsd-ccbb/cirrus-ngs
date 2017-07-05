#!/bin/bash

mkdir /scratch/workspace

sample_name=$1

workspace=/scratch/workspace/$sample_name
rootdir=/shared/workspace/software
samtools=$rootdir/samtools/samtools-1.1/samtools

mkdir -p /scratch/workspace/$sample_name/logs

if [ ! -f $workspace/$sample_name"_1.fq.gz" ]; then
   aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/merge_fastq/$sample_name/$sample_name"_1.fq.gz" $workspace/
   gunzip $workspace/$sample_name"_1.fq.gz"
fi


if [ ! -f $workspace/$sample_name"_2.fq.gz" ]; then
   aws s3 cp s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/merge_fastq/$sample_name/$sample_name"_2.fq.gz" $workspace/
   gunzip $workspace/$sample_name"_2.fq.gz"
fi

java -jar /shared/workspace/software/Trimmomatic-0.30/trimmomatic-0.30.jar PE -threads 8 -phred33 -trimlog /shared/workspace/software/Trimmomatic-0.30/trimlog.log $workspace/$sample_name"_1.fq" $workspace/$sample_name"_2.fq" $workspace/$sample_name"_1.trim.fq" $workspace/$sample_name"_1.unpaired.fq" $workspace/$sample_name"_2.trim.fq" $workspace/$sample_name"_2.unpaired.fq" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm -r $workspace/$sample_name"_1.fq"
rm -r $workspace/$sample_name"_2.fq"
