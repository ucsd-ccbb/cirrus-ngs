#!/bin/bash

export PATH=/shared/workspace/software/anaconda2/bin/:$PATH

ref_genes=/shared/workspace/software/gencode/gencode.v19.annotation.gtf
samtools=/shared/workspace/software/samtools/samtools-1.1/samtools

input_file=$1
output_file=$2

if [ ! -f $output_file ]; then
   $samtools view $input_file | sort -s -k 1,1 - |  htseq-count - $ref_genes > $output_file
fi
exit 0
