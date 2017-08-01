#!/bin/bash

bowtie=/shared/workspace/software/bowtie-1.0.1/bowtie
bowtie_index=/shared/workspace/software/bowtie_index/hsapiens_hg19/genome

input_fastq=$1

if [ ! -f "$input_fastq.sam" ]; then
   $bowtie $bowtie_index $input_fastq > $input_fastq.sam
fi

exit 0
