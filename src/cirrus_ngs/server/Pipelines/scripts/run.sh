#!/bin/bash

yaml_file=$1
log_dir=$2
pipeline=$3     #specific pipeline file

echo $yaml_file
echo $log_dir
echo $pipeline

mkdir -p $log_dir
log_file=$log_dir/'run.log'
exec 1>>$log_file
exec 2>>$log_file

/shared/workspace/software/anaconda3/bin/python /shared/workspace/Pipelines/$pipeline $yaml_file
