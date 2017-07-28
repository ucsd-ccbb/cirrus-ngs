#!/bin/bash

yaml_file=$1
log_dir=$2

mkdir -p $log_dir

exec 1>>$log_dir/run.log
exec 2>>$log_dir/run.log

/shared/workspace/software/anaconda3/bin/python /shared/workspace/WGSPipeline/WGSPipeline.py $yaml_file

echo "Clearing /scratch/ ..."

list=(`ls /scratch/`)

for file in "${list[@]}";
do
    if [[ $file != "lost+found" ]];
    then
        echo $file
        rm -r /scratch/$file
    fi
done


#rm -r /scratch/`basename $project_name`

echo "Data has been cleared"
echo
