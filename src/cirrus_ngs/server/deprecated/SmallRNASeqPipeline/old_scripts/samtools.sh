#!/bin/sh

/shared/workspace/software/samtools/samtools-1.1/samtools merge -h $headSAMFile $outputBAMFile $inputFiles

/shared/workspace/software/samtools/samtools-1.1/samtools view $fastqFile".bam" > $fastqFile".sam"

rm $fastqFile"."?".bam"

echo "done."
