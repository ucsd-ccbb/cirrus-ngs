#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
normal_sample=$5
tumor_sample=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
num_threads=${11}
chromosome=${12}
##TUMOR VS NORMAL


#logging
log_dir=$log_dir/$normal_sample
mkdir -p $log_dir
log_file=$log_dir/"mutect.$chromosome.log"
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$normal_sample
mkdir -p $workspace
mkdir -p $workspace/temp

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME"_$chromosome" $status_file

##DOWNLOAD##
if [ ! -f $workspace/$normal_sample$file_suffix ] || [ ! -f $workspace/$normal_sample$file_suffix.bai ] \
    || [ ! -f $workspace/$tumor_sample$file_suffix ] || [ ! -f $workspace/$tumor_sample$file_suffix.bai ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    aws s3 cp $input_address/$normal_sample/$normal_sample$download_suffix $workspace/
    aws s3 cp $input_address/$normal_sample/$normal_sample$download_suffix.bai $workspace/
    gunzip -q $workspace/$normal_sample$download_suffix

    aws s3 cp $input_address/$tumor_sample/$tumor_sample$download_suffix $workspace/
    aws s3 cp $input_address/$tumor_sample/$tumor_sample$download_suffix.bai $workspace/
    gunzip -q $workspace/$tumor_sample$download_suffix

    #always download forward reads
fi
##END_DOWNLOAD##

##MUTECT##

if [ -z $cosmic ]
then
    references="--dbsnp $dbsnp"
else
    references="--dbsnp $dbsnp --cosmic $cosmic"
fi

check_exit_status "$java -Djava.io.tmpdir=$workspace/temp -Xmx4g -jar $gatk \
    -T MuTect2 \
    -nct $num_threads \
    -R $genome_fasta \
    -I:normal $workspace/$normal_sample$file_suffix \
    -I:tumor $workspace/$tumor_sample$file_suffix \
     $references \
    -L $chromosome \
    -o $workspace/$normal_sample'_vs_'$tumor_sample.$chromosome.vcf" $JOB_NAME"_$chromosome" $status_file
##END_MUTECT##


##UPLOAD##
out_file=$normal_sample'_vs_'$tumor_sample.$chromosome.vcf
aws s3 cp $workspace $output_address --exclude "*" --include "$out_file" --recursive
##END_UPLOAD##
