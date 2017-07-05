#!/bin/bash

workspace=/scratch/workspace
rootdir=/shared/workspace/software

if [ ! -d $workspace ]; then
    mkdir $workspace
fi

#used in header id
groupID=$1
phenotype = $2


