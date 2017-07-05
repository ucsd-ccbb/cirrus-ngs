#!/bin/bash

workspace=/scratch/workspace

if [ ! -d $workspace ]; then
        mkdir $workspace
fi

fileName=$1
s3URL=$2

# redirecting all output to a file
exec 1>>$3/$HOSTNAME"_s3download.o"
exec 2>>$3/$HOSTNAME"_s3download.o"

if [ ! -f $workspace/$fileName ]; then
   aws s3 cp $s3URL/$fileName $workspace
fi

