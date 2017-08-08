#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
makeTagDirectory=/shared/workspace/software/homer/bin/makeTagDirectory
output_tag_folder=$1
input_file=$2

$makeTagDirectory $output_tag_folder $input_file
