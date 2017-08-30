#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
findPeaks=/shared/workspace/software/homer/bin/findPeaks

chip_tag_folder=$1
style=$2
output_peak_file=$3
input_tag_folder=$4

mkdir -p $output_peak_file
rm -r $output_peak_file

if [ "$control_tag_folder" == " " ] ; then
   #$findPeaks $input_tag_folder -region -size 1000 -minDist 2500 -style $style -o $output_peak_file
   $findPeaks $chip_tag_folder -region -size 1000 -minDist 2500 -o $output_peak_file
else
   $findPeaks $chip_tag_folder -region -size 1000 -minDist 2500 -o $output_peak_file -i $input_tag_folder
fi
