#!/bin/bash

dataFolder=$1
S3UploadURL=$2

echo "uploading file: "$dataFolder

aws s3 cp $dataFolder $S3UploadURL/ --recursive --include "*" --exclude "*.fastq" 
