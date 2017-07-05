#!/bin/bash

inputFile=$1
S3DownloadURL=$2
dataFolder=$3

echo "downloading file: "$S3DownloadURL/$inputFile

if [ ! -f $dataFolder$inputFile ]
then
    aws s3 cp $S3DownloadURL/$inputFile $dataFolder/ 
fi
