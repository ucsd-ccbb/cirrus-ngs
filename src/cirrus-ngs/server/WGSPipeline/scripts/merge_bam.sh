#!/bin/bash

workspace=/scratch/workspace
rootdir=/shared/workspace/software

outputFile=$1

# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_mergebam.o"
exec 2>>$2/$HOSTNAME"_mergebam.o"

sambamba=$rootdir/sambamba/0.4.7/bin/sambamba
numThreads=16

# store arguments in a special array 
args=("$@")
# get number of elements 
ELEMENTS=${#args[@]}

for (( i=2;i<$ELEMENTS;i++)); do
    bamFile="$workspace/${args[${i}]} "
    fileList=$fileList$bamFile
done

if [ ! -f $workspace/$outputFile.final.bam ]; then

   $sambamba merge -t $numThreads $workspace/$outputFile.final.bam $fileList

fi

