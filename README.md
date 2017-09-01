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

## The Design File
This txt file specifies what samples will be used in this project.
It supports both a two column and three column tab-separated format

#### Two column format
In the two column format the first column is the filename of the sample  
   * If the sample is paired end the first column should be:  
   * `name_of_forward_end_file,name_of_reverse_end_file`
   * Note that the two files are only separated by a comma, no spaces  
The second column is the name of the group associated with that sample  
Group names are used for variant calling. Samples with the same group will have their vcf files merged and the group-based vcf files will be compared to one another.  

#### Three column format
The three column format has the same first two columns as the two column format.  
The third column is an identifier that is either from Normal, Tumor, Chip, or Input (_case sensitive_)  
   * The Normal/Tumor identifiers are used for mutect in the WGS pipeline  
   * The Chip/Input indentifiers are used throughout the ChipSeq pipeline 
  
If two files form a Normal/Tumor or Chip/Input pair they must have the same group and directly follow one another  
Also, each group in the three column format must have exactly one of each identifier (one Normal && one Tumor) || (one Chip && one Input)  

##### Examples:

##### GOOD
```
sample1_forward,sample1_reverse<TAB>group1<TAB>Normal
sample2_forward,sample2_reverse<TAB>group1<TAB>Tumor
```
```
sample1<TAB>group1<TAB>Chip
sample2<TAB>group1<TAB>Input
```

##### BAD
```
sample1<TAB>group1<TAB>Normal
sample2<TAB>group2<TAB>Tumor
```
```
sample1<TAB>group1<TAB>Normal
sample2<TAB>group1<TAB>Input
```
```
sample1<TAB>group1<TAB>Tumor
sample2<TAB>group1<TAB>Tumor
```
```
sample1<TAB>group1<TAB>Tumor
sample2<TAB>group1<TAB>Normal
```
^ wrong order of samples in last one  


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
        * **.fq or .fastq if uses_chromosomes is false**
        * .trim{} -> .trim(.fq|.fastq)
        * **${curr_chromosome_number} if uses_chromsomes is true**
        * .final.{}.bam -> .final.$CHROMOSOMENUMBER.bam
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

By default the tool will be run on all the samples in the project. Each tool can also be run on all samples at once, on each group, and on pairs of samples if needed.

#### Pipeline specific yaml files
