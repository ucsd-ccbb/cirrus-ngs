#!/bin/bash

bowtie=/shared/workspace/software/bowtie-1.0.1/bowtie
bowtie_index=/shared/workspace/software/bowtie_index/hairpin_human/hairpin_human

>&2 echo "bowtie output file: " $fastqFile".sam"

if [ ! -f $fastqFile".sam" ]; then
   $bowtie $bowtie_index -n 1 -l 8 -S $fastqFile > $fastqFile".sam"
fi

exit 0
