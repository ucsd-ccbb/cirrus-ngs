#!/bin/bash

file_name=$1
workspace=/mnt/efs/workspace/$file_name
rootdir=/shared/workspace/software

mkdir -p $workspace/logs

# redirecting all output to a file
exec 1>>$workspace/logs/"bam2fastq.o"
exec 2>>$workspace/logs/"bam2fastq.e"

numThreads=32
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba

tmpDir=$workspace/temp

echo "Sort is processing: $workspace/$file_name.bam"

### Sort BAM File ####
if [ ! -f $workspace/$file_name.sort.bam ]; then
   $sambamba sort -t $numThreads -m 4G -n --tmpdir $tmpDir -o $workspace/$file_name.sort.bam $workspace/$file_name.bam
fi

echo "converting bam to fastq...."
if [ ! -f $workspace/$file_name'_R1.fq' ] && [ ! -f $workspace/$file_name'_R2.fq' ]; then
   /shared/workspace/software/bedtools2/bin/bamToFastq -i $workspace/$file_name.sort.bam -fq $workspace/$file_name'_R1.fq' -fq2 $workspace/$file_name'_R2.fq'

   echo "renaming the raw bam file."
   mv $workspace/$file_name.bam $workspace/$file_name.raw.bam	
   mv $workspace/$file_name.sort.bam $workspace/$file_name.raw.sort.bam
fi

