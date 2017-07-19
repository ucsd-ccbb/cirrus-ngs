#!/bin/sh

yamlFile=$1

java -jar /shared/workspace/SmallRNASeqPipeline/scripts/SmallRNASeqPipeline.jar $yamlFile > /shared/workspace/SmallRNASeqPipeline/nohup.out

