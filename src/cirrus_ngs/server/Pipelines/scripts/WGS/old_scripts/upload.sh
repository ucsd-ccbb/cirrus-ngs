#!/bin/bash

sample_file=$1

aws s3 cp /scratch/workspace/$sample_file/ s3://analysis-data-by-ccbb/20170320_Shuling_Califano_dnaseq/analysis_results/$sample_file/ --recursive

rm -r /scratch/workspace/$sample_file/

