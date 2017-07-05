#!/bin/bash

workspace=/scratch/workspace

fileName=$1
# redirecting all output to a file
exec 1>>$2/$HOSTNAME"_compressfastq.o"
exec 2>>$2/$HOSTNAME"_compressfastq.o"

if [ ! -f $workspace/$fileName.gz ]; then
   gzip $workspace/$fileName
fi

