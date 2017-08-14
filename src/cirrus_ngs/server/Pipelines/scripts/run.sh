#!/bin/bash

yaml_file=$1
log_dir=$2
pipeline=$3     #specific pipeline file

mkdir -p $log_dir

exec 1>>$log_dir/run.log
exec 2>>$log_dir/run.log

/shared/workspace/software/anaconda3/bin/python /shared/workspace/Pipelines/$pipeline $yaml_file
