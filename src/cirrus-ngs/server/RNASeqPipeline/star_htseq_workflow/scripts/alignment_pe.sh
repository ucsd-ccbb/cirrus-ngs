#!/bin/bash

samtools=/shared/workspace/software/samtools/samtools-1.1/samtools
STAR=/shared/workspace/software/STAR/2.5.1a/STAR/bin/Linux_x86_64/STAR
genomeDir=/shared/workspace/software/genomes/Hsapiens/star

input_fastq=$1
output_file=$2

mkdir -p $output_file

if [ ! -f $output_file/"Aligned.out.sam" ]; then
   $STAR --genomeDir $genomeDir --readFilesIn $input_fastq"_1.trim.fastq" $input_fastq"_2.trim.fastq" --outFileNamePrefix $output_file/
fi

if [ ! -f $output_file/"Aligned.out.bam" ]; then
   $samtools view -Sb $output_file/"Aligned.out.sam" > $output_file/"Aligned.out.bam"
fi

if [ ! -f $output_file/"Aligned.out.sorted.bam" ]; then
   $samtools sort -m 2G -@ 4 $output_file/"Aligned.out.bam" $output_file/"Aligned.out.sorted"
fi

if [ ! -f $output_file/"Aligned.out.sorted.bam.bai" ]; then
   $samtools index $output_file/"Aligned.out.sorted.bam"
fi
exit 0
