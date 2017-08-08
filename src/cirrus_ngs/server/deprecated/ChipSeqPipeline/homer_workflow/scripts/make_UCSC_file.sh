#!/bin/bash

export PATH=$PATH:/shared/workspace/software/homer/bin
makeUCSCfile=/shared/workspace/software/homer/bin/makeUCSCfile
input_tag_folder=$1

$makeUCSCfile $input_tag_folder -o auto
