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

#logging
log_dir=$log_dir/$group_name
mkdir -p $log_dir
log_file=$log_dir/'combine_vcf.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$group_name
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

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

    #download all separated vcf and bam files
    for file in $files_in_group
    do
        aws s3 cp $input_address/$file/$file$file_suffix $workspace/
    done
fi
##END_DOWNLOAD##

variant_list=""

for file in $files_in_group
do
    variant_list=$variant_list"--variant $workspace/$file$file_suffix "
done


##COMBINEVCF##
check_exit_status "$java -Xmx2g -Djava.io.tmpdir=$workspace/temp \
    -jar $gatk \
    -T CombineGVCFs \
    -R $genome_fasta \
    $variant_list \
    -o $workspace/$group_name.merged.g.vcf" $JOB_NAME $status_file

check_exit_status "$java -Xms454m -Xmx3181m -Djava.io.tmpdir=$workspace/temp \
    -jar $gatk \
    -T GenotypeGVCFs \
    -A QualByDepth -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio -A DepthPerSampleHC -A InbreedingCoeff \
    -R $genome_fasta \
    -nt $num_threads \
    --variant $workspace/$group_name.merged.g.vcf \
    -o $workspace/$group_name.vcf \
    --dbsnp $dbsnp" $JOB_NAME $status_file

check_exit_status "check_outputs_exist $workspace/$group_name.vcf" $JOB_NAME $status_file

#END_COMBINEVCF##

##UPLOAD##
aws s3 cp $workspace/ $output_address --exclude "*" --include "$group_name.vcf" --recursive
##END_UPLOAD
