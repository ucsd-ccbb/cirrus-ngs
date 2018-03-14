# cirrus-ngs
Cloud-optimized next generation sequencing primary analysis pipelines for whole genome and exome variant calling, RNA-seq, miRNA-seq, ChIP-seq, DNA-seq.

## Dependencies
All dependencies can be installed with pip
* paramiko
* cfncluster
* yaml
* scp
* aws-cli

## Pipeline and Workflow - Explained
As a convention, we say "pipeline" is at a higher level than "workflow". 
Pipeline corresponds to a general type of sequencing, such as WGS Pipeline, which stands for Whole Genome Sequencing Pipeline. 
Meanwhile, each pipeline is allowed to have multiple "workflows". 
For instance, RNA-Seq Pipeline contains four different workflows (star_gatk, star_htseq, star_rsem, kallisto), each with its own steps.


## Supported Pipelines
* WGSPipeline 
* RNASeqPipeline
* ChiPSeqPipeline
* SmallRNASeqPipeline





## General Overview
First the user creates a design file (format described [below](#design)). The jupyter notebook for the user's chosen pipeline requires such a design file and multiple parameters specified within the first cell. The notebook creates a yaml file summarizing all of the user input and transfers that file to the cluster. Cluster-native code then uses that yaml file, along with multiple configuration files (described [below](#config)), to sequentially execute the analysis steps specified by the user in a distributed fashion. Upon completion of each step the output will be uploaded to the user's s3 output bucket and can be accessed at any point.    

#### A rough diagram:
```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||||||||@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Local                                                         transfer                                           Remote
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||||||||@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
DesignFile + Parameters -> Notebook -> PipelineManager -> yaml ------> Pipeline + config -> shell scripts foreach step  
                                                                                            /\                      ||
                                                                                            ||                      || 
                                                                                            ||Download        Upload||
                                                                                            ||precursors     outputs||
                                                                                            ||                      \/
                                                                                          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                                                                          @          USER S3          @
                                                                                          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
```

#### Output S3 Address Organization 
One parameter the user supplies is the s3_output_files_address, which is used as a base directory for all project output.  
Output files can be found at `$s3_output_files_address/$project_name/$workflow/$sample_name` where $sample_name is the name of the forward read file. Other variables are specified by the user in the notebook. All output for a sample will be under the directory with the name of its forward read file.   
__Note__: After alignment all output files associated will a sample will take on the name of the forward read file.   
&nbsp;&nbsp;Example: If Sample is made up of Sample_R1.fq and Sample_R2.fq alignment will output Sample_R1.sam or Sample_R1.bam.
    
## The Design File <a name="design"></a>
This txt file specifies what samples will be used in this project.
It supports both a two column and three column tab-separated format. 

### Two column format
In the two column format, the __first column__ is the filename of the sample (with extensions: e.g. fastq.gz). 
   * If the sample is paired end, the first column should be:  
   * `forward_end_file,reverse_end_file`
   * Note that the two files are only separated by a comma, __no spaces__  
   
The __second column__ is the name of the group associated with that sample. 
Group names are up to the user, but they are generally set to experimental conditions. For example, all control samples
can be given a group named "control". 


#### Examples
For paired-end:
```
sample1_forward.fastq.gz,sample1_reverse.fastq.gz<TAB>groupA
sample2_forward.fastq.gz,sample2_reverse.fastq.gz<TAB>groupB
```
For single-end:
```
sample1.fastq.gz<TAB>groupA
```

### Three column format
The three column format has the same first two columns as the two column format.  
The third column is an identifier that is either Normal, Tumor, Chip, or Input (_case sensitive_)  
   * The Normal/Tumor identifiers are used for somatic variant calling in the WGS pipeline  
   * The Chip/Input indentifiers are used throughout the ChiPSeq pipeline (Input is used to normalize)
  
If two files form a Normal/Tumor or Chip/Input pair they must have the same group.
Also, each group in the three column format must have exactly one of each identifier (one Normal && one Tumor) || (one Chip && one Input). Do not mix Chip/Input and Normal/Tumor pairs in the same design file.

#### Examples:

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


## Adding additional tools
All tools are run by Pipeline.py on the cluster. Adding additional tools requires:
1. A shell script following the standard format
2. An entry in the workflow specific yaml file

## The shell script format 
```bash
#!/bin/bash

project_name=$1
workflow=$2
file_suffix=$3      #extension of input file, does not include .gz if present in input
root_dir=$4
fastq_end1=$5       #forward reads
fastq_end2=$6       #reverse reads (can be NULL)
input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
log_dir=$9          #directory for log files
is_zipped=${10}     #either "True" or "False", indicates whether input is gzipped
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
workspace=$root_dir/$project_name/$workflow/$fastq_end1
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
Shell scripts used to call tools follow a specific format. The first 10 arguments exist in every shell script. Additional arguments can be specified through the configuration files. 

## Configuration Files <a name="config"></a>
#### Workflow specific yaml files
Each workflow has its own configuration file with two specific sections.  
First, there must be a steps list that contains the order of the steps that can be run. __Order is enforced__. Whatever
order you lists the steps in is what order they will be called. Don't list steps in the wrong order (e.g. haplotyping before
aligning), because the steps will be called in whatever order is listed.

  Example: 
  ```yaml
  steps:
    - "fastqc"
    - "trim"
    - "bwa"
    - "sort"
    - "dedup"
    - "split"
    - "postalignment"
    - "haplotype"
    - "mutect"
    - "merge"
    - "bam_merge"
    - "pair_vcf_merge"
    - "combine_vcf"
    - "filter"
  ```
  
Second, there is a series of configuration entries for each step in the given workflow with the following format.
```yaml
NAME_OF_TOOL:
    script_path: "NAME_OF_SHELL_SCRIPT"
    download_suffix: "NAME_OF_FILE_EXT"
    input_is_output: BOOL_VALUE
    can_be_zipped: BOOL_VALUE
    uses_chromosomes: BOOL_VALUE
    extra_bash_args:
      - LIST
      - EXTRA
      - ARGS
      - HERE
```
NAME_OF_TOOL is what the user will input into analysis_steps in the jupyter notebook in order to call this tool.

Basic notes about the required entries for each tool:
* script_path
    * the name of the shell script to be executed
    * should not contain any file extensions
    * can contain a path to the shell script relative to the /shared/workspace/Pipelines/scripts/ directory
* download_suffix:
    * file extension of prerequisite files for this step
    * if ~ is passed in the prereq extension will default to the extension of the samples for this project (.fq or .fastq)
    * this should contain the "." in the extension (".fq" not "fq")
    * can contain one {} format specifier that takes on a value of
        * **.fq or .fastq if uses_chromosomes is false**
        * .trim{} -> .trim(.fq|.fastq)
        * **${curr_chromosome_number} if uses_chromsomes is true**
        * .final.{}.bam -> .final.$CHROMOSOMENUMBER.bam
* input_is_output
    * boolean value describing location of prerequisite files for this step
    * if true prerequisites will be downloaded from user's specified output s3 bucket instead of their input s3 bucket
    * most steps pass in True, only steps that operate solely on fastq files pass in False
* can_be_zipped
    * boolean value describing if the prerequisites for this step can exist in gzipped format
    * when set to true the download will be for "download_suffix.gz" if the original samples are gzipped
    * when set to false the .gz extension will never be set
    * Note: this is only relevant when a tool is gzipped and the step needs to operate on an unzipped version
      * If a step needs .fastq files and the input has a .fastq.gz, pass in .fastq to download_suffix and True to can_be_zipped
      * If a step needs .vcf.gz files (and doesn't want to unzip them) pass in .vcf.gz to download_suffix and False to can_be_zipped
* uses_chromosomes
    * boolean value describing if this step should be run on each chromosome
      * Precursors must have been split by chromosome
    * the last bash argument to the script will be the current chromosome number
    * shell script arguments:
    ```bash
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
    EXTRA_BASH_ARGUMENTS
    chromosome=${last_parameter}
    ```
* extra_bash_args
  * a list of extra arguments to pass into the shell script for this step
  * can be an empty list if no extra arguments are needed (`extra_bash_args: []`)
  * **If extra parameters are required the first one must be a integer number of threads needed for that step**
  * shell script for the step must support these extra arguments

By default the tool will be run thon all the samples one at a time in the project. Each tool can also be run on all samples at once, on each group, and on pairs of samples if needed. An additional field can be added to each tool's specification to force a different kind of run. **Only one such extra field should be added to each step**.

* all_samples:
  * If set to true the tool will be run once on all samples within the project. These samples will be passed into the tool's shell script as a space-delimited list of the sample forward-read file names. 
  * The input and output addresses will not contain the sample name as the others do. Intead of `$path/proj_name/workflow/sample` the output will be to `$path/proj_name/workflow`
  * shell script arguments:
  ```bash
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
  all_samples=${11}  #this isn't an extra bash arg, it's created for you by the pipeline
  EXTRA_BASH_ARGUMENTS
  ```
  
* by_pair:
  * If set to true the tool will be run on pairs of samples. Sample pairs are determined by the user's design file. Pair-based analysis requires the third field in the design file. Samples that constitute a pair must be the only two samples in their group. More details in design file section.
  * The output address parameter will be set to `$path/proj_name/workflow/normal_sample_name`, where the normal sample was the sample specified as either Chip or Normal. 
  * The input address parameter is partially determined by input_is_output; if input_is_output is true then input address is set to the `$user_output_address/proj_name/workflow`. This allows for more control in downloading files from different s3 buckets. If input_is_output is false then the input address will be the user's specified input address.
  * If by_pair is set to true for some tool but the design file doesn't have a third field then said tool will be run on a by sample basis instead.  
    * The argument generation will be on a by sample basis, but the shell script must support such a change.
  * shell script arguments:
  ```bash
  project_name=$1
  workflow=$2
  file_suffix=$3  #extension of input file, does not include .gz if present in input
  root_dir=$4
  first_in_pair=$5
  second_in_pair=$6
  input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
  output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
  log_dir=$9
  is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
  EXTRA_BASH_ARGUMENTS
  ```
  
* by_group:
  * If set to true the tool will be run on each group of samples as specified by the design file. Group based analysis will create a new directory `$path/$proj_name/$workflow/$group_name` to store analysis for that group. 
  * When run by_group the shell script will have an argument containing a space-delimited list of samples in that group.
  * The output address and input address parameters will be set in the same manner as the by_pair output and input addresses. However, instead of `$path/proj_name/workflow/normal_sample_name` the output address will be set to `$path/proj_name/workflow/group_name` 
  * shell script arguments:
  ```bash
  project_name=$1
  workflow=$2
  file_suffix=$3  #extension of input file, does not include .gz if present in input
  root_dir=$4
  group_name=$5
  fastq_end2=$6     #this is always "NULL" for by_group iteration
  input_address=$7    #this is an s3 address e.g. s3://path/to/input/directory
  output_address=$8   #this is an s3 address e.g. s3://path/to/output/directory
  log_dir=$9
  is_zipped=${10}    #either "True" or "False", indicates whether input is gzipped
  group_samples=${11} #not an "extra" argument, this is generated by the pipeline for you
  EXTRA_BASH_ARGUMENTS
  ```



#### software.conf
This configuration file contains environment variables for paths to executables and reference files. If more variables are required they can simply be added to this file anywhere. When writing extra shell scripts in cirrus one shouldn't rely on any
relative paths. The computing nodes do not inherit the root directory of the head node. Consider installing samtools
and soft linking /usr/bin/samtools to the samtools executable on our head node. None of the computing nodes will be able to
find "samtools" because they will not inherit the /usr/bin/samtools soft link. Instead, samtools is installed under the /shared
directory and an environment variable "$samtools" is defined in software.conf that contains the absolute path to the samtools executable. This way both head and computing nodes can use the tool with $samtools.

This file also contains three bash functions used in all shell scripts to prevent duplicate tool running.
* check_step_already_done
  * args: name of the job being run, path to status.log 
  * This function checks the status.log file passed in to see if the currently running job has already passed. If it has already passed, then the shell script will terminate early and not run the step unnecessarily.

* check_exit_status
  * args: string command to be run, name of the job being run, path to status.log
  * This function checks the exit status of the command passed in. If it is non-zero it attempts the step twice more. If it is still non-zero then the error messages are logged and a "failed" specifier is added to status.log
  
* check_outputs_exist
  * args: all outputs to check existence of (outputs can be directories or files)
  * returns 1 if any of the given outputs do not exist, 0 otherwise
  

## Cluster organization
Clusters are created using an EBS snapshot containing all the tools (bowtie, samtools, etc) and software needed to run cirrus-ngs. These clusters follow a very specific directory organization which must be maintained if the user wants to add anything to the cluster. All relevant files and directories are rooted at /shared/workspace/ which will be shortened to $root in the following sections. 
### Directories
* $root/logs/
* $root/software/
* $root/Pipelines/

#### $root/logs/
This directory contains the logs from any pipeline. A given pipeline's logs will be in $root/logs/$pipeline/$workflow/$project_name. The $log_dir parameter in each shell script contains this path. The user doesn't need to set it, it's handled already.

#### $root/software/
This directory contains all the tools that cirrus relies on to run. Each tool is under its own directory.

**Example:**
```$root/software/$tool_name/$version/```

A tool can have multiple versions. Under the version is all the software associated with any given tool. For example, currently
the STAR tool has 3 versions: 2.3.0e, 2.5.1a, and 2.5.3a. The STAR directory is organized as follows:

* $root/software/STAR/
  * 2.3.0e/
    * <all STAR 2.3.0e files>
  * 2.5.1a/
    * <all STAR 2.5.1a files>
  * 2.5.3a/
    * <all STAR 2.5.3a files>
    
The software directory also contains a references directory ($root/software/references)
##### $root/software/references/

This directory contains all of the reference files needed for cirrus. This includes, but isn't limited to, fasta references, gtf files, and alignment indices (for different alignment tools). The directory is organized as follows:

* $root/software/references/
  * $organism_name/ (Hsapiens, Mmusculus, etc.)
    * $assembly_name/ (hg19, etc.)
      * annotations/
        * gtf files, general gene annotations
      * indices/
        * $alignment_tool_name/ (STAR, bowtie, etc.)
          * index built for that alignment tool
      * sequence/
        * fasta files for sequence references
      * variation/
        * vcf references (dbsnp, comsic, etc.)
        
This structure isn't set in stone, but provides an organized way of storing any references thay may be needed. Following is an example of how the references directory would look if we had only the human hg19 references.

* $root/software/references/
  * Hsapiens/
    * hg19/
      * annotations/
        * gencode.v19.annotation.gtf
      * indices/
        * bwa/
          * ucsc.hg19.fasta.amb
          * ucsc.hg19.fasta.ann
          * ucsc.hg19.fasta.bwt
          * ucsc.hg19.fasta.pac
          * ucsc.hg19.fasta.sa
        * STAR/
          * chrLength.txt
          * chrNameLength.txt
          * chrName.txt
          * chrStart.txt
          * exonGeTrInfo.tab
          * exonInfo.tab
          * geneInfo.tab
          * Genome
          * genomeParameters.txt
          * human.chrlist
          * human.grp
          * human.idx.fa
          * humanLog.out
          * human.n2g.idx.fa
          * human.seq
          * human.ti
          * human.transcripts.fa
          * SA
          * SAindex
          * sjdbInfo.txt
          * sjdbList.fromGTF.out.tab
          * sjdbList.out.tab
          * transcriptInfo.tab
        * bowtie/
          * genome.1.ebwt
          * genome.2.ebwt
          * genome.3.ebwt
          * genome.4.ebwt
          * genome.fa
          * genome.rev.1.ebwt
          * genome.rev.2.ebwt
        * bowtie2/
        * kallisto/
          * kallisto_index
      * sequence/
        * ucsc.hg19.dict
        * ucsc.hg19.fasta
        * ucsc.hg19.fasta.fai
      * variation/
        * 1000G_omni2.5.hg19.sites.vcf
        * 1000G_omni2.5.hg19.sites.vcf.idx
        * 1000G_phase1.indels.hg19.sites.vcf
        * 1000G_phase1.indels.hg19.sites.vcf.idx
        * 1000G_phase1.snps.high_confidence.hg19.sites.vcf
        * 1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx
        * cosmic_hg19.vcf
        * dbsnp_138.hg19.vcf
        * dbsnp_138.hg19.vcf.idx
        * hapmap_3.3.hg19.sites.vcf
        * hapmap_3.3.hg19.sites.vcf.idx
        * Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
        * Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx


### $root/Pipelines/

This directory contains the core of the non-local cirrus-ngs implementation. The main file, Pipeline.py, handles generating shell script arguments and submitting jobs based on the user's specifications. The directory is organized as follows. The config directory contains configuration files and information. The scripts directory contains all the shell scripts associated with all possible steps. The yaml_files directory contains the yaml files created from the user's input into the jupyter notebook. The util directory contains auxillary scripts and tools written specifically for cirrus.

* $root/Pipelines/
  * Pipeline.py
  * config/
    * software.conf
    * $pipeline/
      * configuration files for all workflows in $pipeline
      * files are named ${pipeline}_${workflow}.yaml
  * scripts/
    * run.sh (this is the main entry point, it calls Pipeline.py)
    * shared shell scripts (scripts called by every pipeline/workflow)
    * $pipeline/
      * shell scripts called by all workflows in this pipeline
      * $workflow/
        * shell scripts called by this workflow
  * util/
    * auxillary scripts for cirrus
  * yaml_files/
    * $pipeline/
      * $workflow/
        * $project_name.yaml (created from user's input to juptyer notebook)
        
The script_path filed mentioned in the configuration file section contains script paths relative to $root/Pipelines/scripts/. 
