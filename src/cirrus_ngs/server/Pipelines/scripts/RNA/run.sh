#!/bin/bash

yaml_file=$1
python /shared/workspace/RNASeqPipeline/RNASeqPipeline.py $yaml_file > /shared/workspace/RNASeqPipeline/nohup.out &
 
