#!/bin/bash

workspace=/scratch/workspace

fileName=$1
s3URL=$2
# redirecting all output to a file
exec 1>>$3/$HOSTNAME"_s3upload.o"
exec 2>>$3/$HOSTNAME"_s3upload.o"

echo "fileName: $fileName"
echo "s3URL: $s3URL"

if [ -f $workspace/$fileName ]; then
  
   aws s3 cp $workspace/$fileName $s3URL/

   rm $workspace/$fileName
fi

