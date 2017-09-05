#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
chip_sample=$4
input_sample=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
num_threads=${10}

#logging
log_dir=$log_dir/$chip_sample
mkdir -p $log_dir
log_file=$log_dir/'find_motifs_genome.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$chip_sample
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

pair_base_name=$chip_sample

##DOWNLOAD##
if [ ! -f $workspace/$chip_sample$style_ext ] || [ "$pairs_exist" == "True" ]
then
    if [ "$pairs_exist" == "True" ]
    then
        pair_base_name=$chip_sample'_vs_'$input_sample
        input_address=$input_address/$project_name/$chip_sample
    fi

    aws s3 cp $input_address/$pair_base_name$style_ext $workspace/
fi
##END_DOWNLOAD##


##ANNOTATEPEAKS##
export PATH=$PATH:/shared/workspace/software/blat
export PATH=$PATH:/shared/workspace/software/weblogo
export PATH=$PATH:/shared/workspace/software/ghostscript-9.19-linux-x86_64
export PATH=$PATH:/shared/workspace/software/samtools/samtools-1.1
check_exit_status "$find_motifs_genome $workspace/$pair_base_name$style_ext \
    $genome $workspace/motifs_$pair_base_name -size 200 -mask -p $num_threads" $JOB_NAME $status_file
##END_ANNOTATEPEAKS##


##UPLOAD##
aws s3 cp $workspace/motifs_$pair_base_name $output_address/motifs_$pair_base_name --recursive
##END_UPLOAD##
