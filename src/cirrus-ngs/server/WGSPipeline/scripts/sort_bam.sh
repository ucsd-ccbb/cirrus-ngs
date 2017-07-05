#!/bin/bash

file_name=$1

workspace=/scratch/workspace/$file_name
rootdir=/shared/workspace/software
numThreads=8
sambamba=$rootdir/sambamba/0.4.7/bin/sambamba

# redirecting all output to a file
exec 1>>$workspace/logs/"sort_bam.o"
exec 2>>$workspace/logs/"sort_bam.e"

tmpDir=$workspace/temp

echo "Sort is processing: $workspace/$file_name.bam"

### Sort BAM File ####
if [ ! -f $workspace/$file_name.sort.bam ]; then
        $sambamba sort -t $numThreads -m 5G --tmpdir $tmpDir -o $workspace/$file_name.sort.bam $workspace/$file_name.bam
fi

if [ ! -f $workspace/$file_name.sort.bam.bai ]; then
        $sambamba index -t $numThreads $workspace/$file_name.sort.bam $workspace/$file_name.sort.bam.bai
fi


