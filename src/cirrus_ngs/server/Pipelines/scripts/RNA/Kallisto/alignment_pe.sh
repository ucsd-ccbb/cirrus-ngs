#!/bin/bash

input_fastq=$1
output_file=$2

/shared/workspace/software/kallisto_linux-v0.42.1/kallisto quant -i /shared/workspace/software/kallisto_linux-v0.42.1/kallisto_index/kallisto_index -o $output_file -l 50 $input_fastq"_1.trim.fastq" $input_fastq"_2.trim.fastq" 

