#!/bin/bash

inputFile=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_sortbam.o"
exec 2>>$2/$HOSTNAME"_sortbam.o"

workspace=/scratch/workspace
rootdir=/shared/workspace/software
numThreads=16
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba

tmpDir=/scratch/workspace/tmp

echo "Sort shell is processing: $workspace/$inputFile.bam"

### Sort BAM File ####
if [ ! -f $workspace/$inputFile.sort.bam ]; then
	$sambamba sort -t $numThreads -m 4G --tmpdir $tmpDir -o $workspace/$inputFile.sort.bam $workspace/$inputFile.bam
	
	rm $workspace/$inputFile.bam
fi

if [ ! -f $workspace/$inputFile.sort.bam.bai ]; then
	$sambamba index -t $numThreads $workspace/$inputFile.sort.bam $workspace/$inputFile.sort.bam.bai
fi

