#!/bin/bash

input=$1
output=$2

if [ ! -d $output ]; then
    mkdir $output
fi

/shared/workspace/software/FastQC/fastqc --noextract -o $output $input
