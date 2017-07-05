#!/bin/bash

rootdir=/shared/workspace/software
workspace=/shared/workspace

inputBAM=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_bam2fastq.o"
exec 2>>$2/$HOSTNAME"_bam2fastq.o"

if [ ! -f $workspace/$inputBAM'_R1.fq' ] && [ ! -f $workspace/$inputBAM'_R2.fq' ]; then
	java -jar $rootdir/picard-1.96/SamToFastq.jar I=$workspace/$inputBAM.bam F=$workspace/$inputBAM'_R1.fq' F2=$workspace/$inputBAM'_R2.fq'
fi

