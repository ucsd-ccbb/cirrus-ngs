#!/bin/bash

yaml_file=$1
log_dir=$2

mkdir -p $log_dir

touch $log_dir/run.log

exec 1>>$log_dir/run.log
exec 2>>$log_dir/run.log

python /shared/workspace/WGSPipeline/WGSPipeline.py $yaml_file

project_name=${yaml_file//_group*/}

echo "Clearing /scratch/ ..."

rm -r /scratch/`basename $project_name`

echo "Data has been cleared"
