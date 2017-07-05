#!/bin/bash

S3UploadURL=$1
dataFolder=$2

echo "uploading file: "$dataFolder

aws s3 cp $dataFolder $S3UploadURL/ --recursive --include "*" --exclude "*.fastq" 
