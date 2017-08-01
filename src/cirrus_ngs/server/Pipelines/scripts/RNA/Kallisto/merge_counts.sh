#!/bin/bash

workflow=$1
project_name=$2
sample_str=$3

python /shared/workspace/RNASeqPipeline/$workflow/scripts/MergeCountFile.py $workflow $project_name $sample_str
