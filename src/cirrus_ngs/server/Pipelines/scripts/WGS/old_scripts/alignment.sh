#!/bin/bash

file_name=$1
workspace=/scratch/workspace/$file_name
rootdir=/shared/workspace/software

if [ ! -d $workspace ]; then
   mkdir $workspace
fi

mkdir $workspace/logs

# redirecting all output to a file
exec 1>>$workspace/logs/"alignment.o"
exec 2>>$workspace/logs/"alignment.e"

groupID="2017_03_24"
phenotype="hpv"

numThreads=8
bwa=$rootdir/bwa/bwa-0.7.12/bwa
samblaster=$rootdir/samblaster/samblaster
samtools=$rootdir/samtools/samtools-1.1/samtools
genomeSeq=/shared/shuling/resource/bwa/human_g1k_v37.fasta

##### To align paired end fastq files with BWA #####
if [[ -f $workspace/$file_name'_1.trim.fq' && -f $workspace/$file_name'_2.trim.fq' ]]; then
   if [ ! -f $workspace/$file_name.bam ]; then
      $bwa mem -M -t $numThreads -R '@RG\tID:1\tPL:illumina\tPU:'$groupID'\tSM:'$phenotype -v 1 $genomeSeq $workspace/$file_name'_1.trim.fq' $workspace/$file_name'_2.trim.fq' | $samblaster | $samtools view -Sb - > $workspace/$file_name.bam
   fi
fi

rm $workspace/$file_name'_1.trim.fq'
rm $workspace/$file_name'_2.trim.fq'
