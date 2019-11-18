#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3  #extension of input file, does not include .gz if present in input
root_dir=$4
fastq_end1=$5
fastq_end2=$6
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9
is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
num_threads=${11}   # number of threads

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'cal_expression.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$workflow/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

#this is the suffix of the input from s3
download_suffix=$file_suffix

if [ "$is_zipped" == "True" ]
then
    download_suffix=$file_suffix".gz"
fi

##DOWNLOAD##
if [ ! -f $workspace/"Aligned"$file_suffix ] || [ ! -f $workspace/"Aligned"$file_suffix ] 
then
    #always download forward reads
    aws s3 cp $input_address/"Aligned"$download_suffix $workspace/ --quiet
fi
##END_DOWNLOAD##


# RSEM calculate expression
if [ "$fastq_end2" == "NULL" ]
then
    check_exit_status "$rsem --bam -p $num_threads $workspace/'Aligned'$download_suffix \
	$rsem_index $workspace/$fastq_end1" $JOB_NAME $status_file
else
    # paired end
    check_exit_status "$rsem --paired-end --bam -p $num_threads $workspace/'Aligned'$download_suffix \
	$rsem_index $workspace/$fastq_end1" $JOB_NAME $status_file
fi

# RSEM converts to genome bam and sort.
check_exit_status "$rsem_tbam2gbam $rsem_index $workspace/$fastq_end1'.transcript.bam' $workspace/$fastq_end1'.genome.bam'" $JOB_NAME $status_file   
check_exit_status "$sambamba sort --tmpdir=$workspace/tmp -m 2G -t $num_threads $workspace/$fastq_end1'.genome.bam' -o $workspace/$fastq_end1'.genome.sorted.bam'" $JOB_NAME $status_file
check_exit_status "$samtools index $workspace/$fastq_end1'.genome.sorted.bam'" $JOB_NAME $status_file

# perform RSeQC
if [ ! -z "$rRNA_bed" ]
then
    check_exit_status "$rseqc_split_bam -i $workspace/$fastq_end1'.genome.sorted.bam' -r $rRNA_bed -o $workspace/$fastq_end1 > $workspace/$fastq_end1.ribosomal.txt" $JOB_NAME $status_file
    check_exit_status "$rseqc_geneBody_coverage -r $HouseKeepingGenes_bed -i $workspace/$fastq_end1'.genome.sorted.bam' -o $workspace/$fastq_end1" $JOB_NAME $status_file
    check_exit_status "$rseqc_read_distribution -i $workspace/$fastq_end1'.genome.sorted.bam' -r $RefSeq_bed > $workspace/$fastq_end1.read.distribution.txt" $JOB_NAME $status_file

# perform samtools stats, for multiqc purposes
    check_exit_status "check_outputs_exist $workspace/$fastq_end1.genes.results \
    $workspace/$fastq_end1.isoforms.results $workspace/$fastq_end1.stat \
    $workspace/$fastq_end1.ribosomal.txt $workspace/$fastq_end1.read.distribution.txt \
    $workspace/$fastq_end1.geneBodyCoverage.txt" $JOB_NAME $status_file
fi

rm $workspace/'Aligned'$download_suffix
rm $workspace/$fastq_end1'.transcript.bam'

##UPLOAD##
aws s3 cp $workspace $output_address --recursive
##END_UPLOAD##

##CLEAN##
#rm -r $workspace
##END_CLEAN##
