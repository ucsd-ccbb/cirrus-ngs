#!/bin/bash

s3URL=$1
fileName=$2

# redirecting all output to a file
exec 1>>$3/$HOSTNAME"_s3delete.o"
exec 2>>$3/$HOSTNAME"_s3delete.o"

echo "fileName: $fileName"
echo "s3URL: $s3URL"
  
aws s3 rm $s3URL/$fileName


