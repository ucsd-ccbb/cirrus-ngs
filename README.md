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
3. An entry in the `Pipeline`.yaml file for extra arguments to shell script
