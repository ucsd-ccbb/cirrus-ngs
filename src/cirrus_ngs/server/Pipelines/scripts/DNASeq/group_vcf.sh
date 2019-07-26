#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
group_name=$5
fastq_end2=$6       #this is always "NULL" for group-based iteration
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
files_in_group=${11}    #all files in current group
num_threads=${12}
chromosome=${13}

#logging
log_dir=$log_dir/$group_name
mkdir -p $log_dir
log_file=$log_dir/"group_vcf.$chromosome.log"
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$group_name"_"$chromosome
mkdir -p $workspace/temp

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME"_$chromosome" $status_file

##DOWNLOAD##
downloads_needed="False"
for file in $files_in_group
do
    if [ ! -f $workspace/$file$file_suffix ]
    then
        downloads_needed="True"
    fi
done

if [ "$downloads_needed" == "True" ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi

    #download all separated vcf and .tbi files
    for file in $files_in_group
    do
        aws s3 cp $input_address/$file/$file.$chromosome$download_suffix $workspace/ --quiet
        aws s3 cp $input_address/$file/$file.$chromosome$download_suffix".tbi" $workspace/ --quiet
        if [ ! -f $workspace/$file.$chromosome$download_suffix".tbi" ]
        then
            $tabix -p vcf $workspace/$file.$chromosome$download_suffix
        fi
    done
fi
##END_DOWNLOAD##

variant_list=""

for file in $files_in_group
do
    variant_list=$variant_list"--variant $workspace/$file.$chromosome$download_suffix "
done

##GenomicsDBImport##
check_exit_status "$java -jar $gatk \
    GenomicsDBImport \
    $variant_list \
    --genomicsdb-workspace-path $workspace/database_"$chromosome"
    --tmp-dir=$workspace/temp \
    -L $chromosome" $JOB_NAME"_$chromosome" $status_file

# change to the directory of workspace because GenomicsDBImport has 
# to use workspace as the current work directory.
cd $workspace

##GenotypeGVCFs##
check_exit_status "$java -jar $gatk \
    GenotypeGVCFs \
    -R $genome_fasta \
    -V gendb://database_"$chromosome" \
    -O $workspace/$group_name.$chromosome.vcf.gz \
    --tmp-dir=$workspace/temp" $JOB_NAME"_$chromosome" $status_file

check_exit_status "check_outputs_exist $workspace/$group_name.$chromosome.vcf.gz" $JOB_NAME"_$chromosome" $status_file
#check_exit_status "$tabix -p vcf $workspace/$group_name.$chromosome.vcf.gz" $JOB_NAME"_$chromosome" $status_file
#END_COMBINEVCF##

##UPLOAD##
aws s3 cp $workspace/ $output_address --exclude "*" --include "$group_name.$chromosome.vcf.gz*" --recursive --quiet
##END_UPLOAD

##CLEAN##
rm -r $workspace
##END_CLEAN##
