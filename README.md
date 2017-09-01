# cirrus-ngs
Cloud-optimized next generation sequencing primary analysis pipelines for whole genome and exome variant calling, RNA-seq, miRNA-seq, ChIP-seq, and ATAC-seq.

## Dependencies
All dependencies can be installed with pip
* paramiko
* cfncluster
* scp
* aws-cli

## Supported Pipelines
* WGSPipeline
* RNASeqPipeline
* ChipSeqPipeline
* SmallRNASeqPipeline

## Adding additional tools
All tools are run by Pipeline.py on the cluster. Adding additional tools requires:
1. A shell script following the standard format
2. An entry in the tools.yaml file
3. An entry in the Pipeline specific yaml file for extra shells script arguments

## The shell script format 
```bash
#!/bin/bash

project_name=$1
file_suffix=$2  #extension of input file, does not include .gz if present in input
root_dir=$3
fastq_end1=$4
fastq_end2=$5
input_address=$6    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$7   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$8
is_zipped=$9    #either "True" or "False", indicates whether input is gzipped
EXTRA ARGUMENTS HERE

#logging
log_dir=$log_dir/$fastq_end1
mkdir -p $log_dir
log_file=$log_dir/'LOGNAMEHERE.log'
exec 1>>$log_file
exec 2>>$log_file

status_file=$log_dir/'status.log'
touch $status_file

#prepare output directories
workspace=$root_dir/$project_name/$fastq_end1
mkdir -p $workspace

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
date
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

check_step_already_done $JOB_NAME $status_file

##DOWNLOAD##
if [ ! -f $workspace/$fastq_end1$file_suffix ]
then
    #this is the suffix of the input from s3
    download_suffix=$file_suffix

    #changes extension if S3 input is zipped
    if [ "$is_zipped" == "True" ]
    then
        download_suffix=$file_suffix".gz"
    fi  

    #always download forward reads
    aws s3 cp $input_address/$fastq_end1$download_suffix $workspace/
    gunzip -q $workspace/$fastq_end1$download_suffix

    #download reverse reads if they exist
    if [ "$fastq_end2" != "NULL" ]
    then
        aws s3 cp $input_address/$fastq_end2$download_suffix $workspace/
        gunzip -q $workspace/$fastq_end2$download_suffix
    fi  
fi
##END_DOWNLOAD##


##TOOLHERE##
check_exit_status "TOOLCALLHERE" $JOB_NAME $status_file
##END_TOOLHERE##


##UPLOAD##
aws s3 cp $workspace $output_address --exclude "*" --include "GLOBTOINCLUDE" --recursive
##END_UPLOAD##
```
Shell scripts used to call tools follow a specific format. The first 9 arguments exist in every shell script. Additional arguments can be specified through the configuration files. 

## Configuration Files
There are two important configuration yaml files for each tool. 
#### tools.yaml
This contains a comprehensive collection of all possible steps for cirrus-ngs.
Each tool has the following format within tools.yaml
```yaml
NAME_OF_TOOL:
    shell_script: "NAME_OF_SHELL_SCRIPT"
    download_suffix: "NAME_OF_FILE_EXT"
    input_is_output: BOOL_VALUE
    can_be_zipped: BOOL_VALUE
    uses_chromosomes: BOOL_VALUE
```
NAME_OF_TOOL is what the user will input into analysis_steps in the jupyter notebook in order to call this tool.

Basic notes about the required entries for each tool:
* shell_script
    * the name of the shell script to be executed
    * should not contain any file extensions
    * can contain a path to the shell script relative to the /shared/workspace/Pipelines/scripts/ directory
* download_suffix:
    * file extension of prerequissite files for this step
    * if ~ is passed in the prereq extension will default to the extension of the samples for this project (.fq or .fastq)
    * this should contain the "." in the extension (".fq" not "fq")
    * can contain one {} format specifier that takes on a value of
        * .fq or .fastq if uses_chromosomes is false **_.trim{} becomes .trim(.fq|.fastq)_**
        * ${curr_chromosome_number} if uses_chromsomes is true 
        **_.final.{}.bam becomes .final.$CHROMOSOMENUMBER.bam)_**
* input_is_output
    * boolean value describing location of prerequisite files for this step
    * if true downloads will be from user's specified output s3 bucket instead of their input s3 bucket
* can_be_zipped
    * boolean value describing if the prerequisites for this step can exist in gzipped format
    * when set to true the download will be for "download_suffix.gz" if the original samples are gzipped
    * when set to false the .gz extension will never be set
* uses_chromosomes
    * boolean value describing if this step should be run on each chromosome
    * the last bash argument to the script will be the chromosome's number

#### Pipeline specific yaml files
