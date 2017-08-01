#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
pos2bed=/shared/workspace/software/homer/bin/pos2bed.pl

input_peak_region_file=$1
output_bed_file=$2

if [ ! -f $output_bed_file ]; then
   $pos2bed $input_peak_region_file > $output_bed_file
fi
