#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
annotatePeaks=/shared/workspace/software/homer/bin/annotatePeaks.pl

input_peak_region_file=$1
annotated_peak_region_file=$2
output_folder=$3
genome=$4

if [ ! -f $annotated_peak_region_file ]; then
   $annotatePeaks $input_peak_region_file $genome -go $output_folder"_GO" -genomeOntology $output_folder"_Ontology" > $annotated_peak_region_file
fi
