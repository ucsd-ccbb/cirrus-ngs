#!/bin/bash

input=$1
output=$2

mkdir $output

/shared/workspace/software/FastQC/fastqc -o $output $input
