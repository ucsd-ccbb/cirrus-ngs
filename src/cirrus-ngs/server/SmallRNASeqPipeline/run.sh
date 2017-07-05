#!/bin/sh

# redirecting all output to a file
exec 1>>"/shared/workspace/SmallRNASeqPipeline/nohup.out"
exec 2>>"/shared/workspace/SmallRNASeqPipeline/nohup.out"

nohup sh /shared/workspace/SmallRNASeqPipeline/scripts/start.sh $1 &

