#!/bin/bash

fileName=$1
chrom=$2
start=$3
end=$4

# redirecting all output to a file
exec 1>>$5/$HOSTNAME"_splitbam.o"
exec 2>>$5/$HOSTNAME"_splitbam.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software

numThreads=1
samtools=$rootdir/samtools/samtools-1.1/samtools
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba

outputBAM=$fileName.$chrom:$start-$end.sort.split.bam

### Split BAM####
if [ ! -f $workspace/$outputBAM ]; then
	$samtools view -b $workspace/$fileName.sort.bam $chrom:$start-$end > $workspace/$outputBAM
fi

if [ ! -f $workspace/$outputBAM.bai ]; then
        $sambamba index -t $numThreads $workspace/$outputBAM $workspace/$outputBAM.bai
fi
