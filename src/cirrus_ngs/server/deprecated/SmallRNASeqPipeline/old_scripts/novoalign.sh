#!/bin/sh

/shared/workspace/software/Novoalign/novocraft/novoalign -d /shared/workspace/software/Novoalign/novoindex/hairpin_musculus.ndx -f $fastqFile -a TGGAATTCTCGGGTGCCAAGG -r All 1 -l 18 -t 45 -h 90 -o FullNW -o SAM > $outputFile".sam"

/shared/workspace/software/samtools/samtools-1.1/samtools view -S $outputFile".sam" -b -o $outputFile".bam"
