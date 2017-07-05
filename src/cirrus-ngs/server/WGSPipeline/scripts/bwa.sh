#!/bin/bash

workspace=/scratch/workspace
rootdir=/shared/workspace/software

if [ ! -d $workspace ]; then
   mkdir $workspace
fi

groupID=$1
phenotype=$2
fastqFile1=$3
fastqFile2=$4
output=$5

# redirecting all output to a file
exec 1>>$6/$HOSTNAME"_bwa.o"
exec 2>>$6/$HOSTNAME"_bwa.o"

numThreads=16
bwa=$rootdir/bwa/bwa-0.7.12/bwa
samblaster=$rootdir/samblaster/samblaster
samtools=$rootdir/samtools/samtools-1.1/samtools
genomeSeq=$rootdir/genomes/Hsapiens/bwa/ucsc.hg19.fasta

##### To align paired end fastq files with BWA #####
if [[ -f $workspace/$fastqFile1 && -f $workspace/$fastqFile2 ]]; then
   if [ ! -f $workspace/$output.bam ]; then
      $bwa mem -M -t $numThreads -R '@RG\tID:1\tPL:illumina\tPU:'$groupID'\tSM:'$phenotype -v 1 $genomeSeq $workspace/$fastqFile1 $workspace/$fastqFile2 | $samblaster | $samtools view -Sb - > $workspace/$output.bam
		
      rm $workspace/$fastqFile1 $workspace/$fastqFile2
   fi
fi

