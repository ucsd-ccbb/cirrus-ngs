#!/bin/bash

yaml_file=$1
log_dir=$2
pipeline=$3     #specific yaml file for a pipeline

mkdir -p $log_dir
log_file=$log_dir/'run.log'
exec 1>>$log_file
exec 2>>$log_file

source /shared/workspace/cirrus-ngs/src/cirrus_ngs/server/Pipelines/config/software.conf

$python3 /shared/workspace/cirrus-ngs/src/cirrus_ngs/server/Pipelines/Pipeline.py $yaml_file $log_dir $pipeline.yaml
