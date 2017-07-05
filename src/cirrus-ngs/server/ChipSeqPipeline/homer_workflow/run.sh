#!/bin/bash

yaml_file=$1
python /shared/workspace/ChiPSeqPipeline/homer_workflow/ChipSeqPipeline.py $yaml_file > /shared/workspace/ChiPSeqPipeline/nohup.out &
 
