#!/bin/bash

yaml_file=$1
log_dir=$2

# print the parameters
echo "In run.sh"
echo "yaml file: "$yaml_file
echo "log_dir: "$log_dir

exec 1>>$log_dir/run.log
exec 2>>$log_dir/run.log


python /shared/workspace/SmallRNASeqPipeline/miRNAPipeline.py $yaml_file
