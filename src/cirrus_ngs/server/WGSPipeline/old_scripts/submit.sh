#!/bin/bash

yaml_files=$1
log_dir=$2

yaml_files_arr=(`tr "," "\n" <<< $yaml_files`)

mkdir -p $log_dir

touch $log_dir/submit.log

exec 1>>$log_dir/submit.log
exec 2>>$log_dir/submit.log

for yaml_file in ${yaml_files_arr[*]}; do
    qsub -o /dev/null -e /dev/null -pe smp 16 /shared/workspace/WGSPipeline/run.sh $yaml_file $log_dir
done
