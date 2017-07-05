#!/bin/bash

if [ ! -d $workspace ]; then
   mkdir $workspace
fi

fastqFile=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_fastqc.o"
exec 2>>$2/$HOSTNAME"_fastqc.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software

fastqc=$rootdir/FastQC/fastqc

##### To run fastQC #####
if [ -f $workspace/$fastqFile.fq.gz ]; then
   $fastqc $workspace/$fastqFile.fq.gz
fi

