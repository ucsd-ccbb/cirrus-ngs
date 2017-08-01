#!/bin/bash

fileName=$1
chrom=$2

workspace=/scratch/workspace/$fileName
rootdir=/shared/workspace/software

numThreads=1
samtools=$rootdir/samtools/samtools-1.1/samtools
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba

outputBAM=$fileName.$chrom.bam

### Split BAM####
if [ ! -f $workspace/$outputBAM ]; then
	$samtools view -b $workspace/$fileName.dedup.bam $chrom > $workspace/$outputBAM
fi

if [ ! -f $workspace/$outputBAM.bai ]; then
        $sambamba index -t $numThreads $workspace/$outputBAM $workspace/$outputBAM.bai
fi


