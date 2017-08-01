#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
export PATH=$PATH:/shared/workspace/software/blat
export PATH=$PATH:/shared/workspace/software/weblogo
export PATH=$PATH:/shared/workspace/software/ghostscript-9.19-linux-x86_64
export PATH=$PATH:/shared/workspace/software/samtools/samtools-1.1

findMotifsGenome=/shared/workspace/software/homer/bin/findMotifsGenome.pl

input_peak_region_file=$1
output_motif_folder=$2
genome=$3

$findMotifsGenome $input_peak_region_file $genome $output_motif_folder -size 200 -mask
