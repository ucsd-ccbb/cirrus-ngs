#!/bin/bash

yaml_file=$1

# redirecting all output to a file
exec 1>>"/shared/workspace/WGSPipeline/nohup.out"
exec 2>>"/shared/workspace/WGSPipeline/nohup.out"

nohup java -jar /shared/workspace/WGSPipeline/WGSPipelineLauncher.jar $yaml_file > /shared/workspace/WGSPipeline/nohup.out &

