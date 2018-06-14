#!/bin/bash

yaml_file=$1
log_dir=$2
pipeline=$3     #specific yaml file for a pipeline

mkdir -p $log_dir
log_file=$log_dir/'run.log'
exec 1>>$log_file
exec 2>>$log_file

source /shared/workspace/Pipelines/config/software.conf

# get all stack information 
aws cloudformation describe-stacks \
    --stack-name $(cut -f2 -d'=' /opt/cfncluster/cfnconfig | head -n1) \
    > stack.json

# extra instance type for compute nodes
instance_type=$(python - <<END
import json

f = open("stack.json")
data = json.load(f)
params = data["Stacks"][0]["Parameters"]

for param in params:
    if param["ParameterKey"] == "ComputeInstanceType":
        print(param["ParameterValue"])
END
)

rm stack.json

# convert instance type to number of cores, save in max_cores for subprocess use
max_cores=$(awk -v instance="$instance_type" -F '[\t]' '$2 == instance {print $5}' $instances_tsv | awk '{print $1}')
export max_cores=$max_cores

$python /shared/workspace/Pipelines/Pipeline.py $yaml_file $log_dir $pipeline.yaml
