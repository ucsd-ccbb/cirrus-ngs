#!/bin/bash

file_name=$1
workspace=/scratch/workspace/$file_name

# redirecting all output to a file
exec 1>>$workspace/logs/"markduplicates.o"
exec 2>>$workspace/logs/"markduplicates.e"

echo "$(date): running picard..."
## Picard
java -jar -Djava.io.tmpdir=$workspace/temp -Xms250m -Xmx20g /shared/workspace/software/picard-1.96/MarkDuplicates.jar INPUT=$workspace/$file_name.sort.bam OUTPUT=$workspace/$file_name.dedup.bam METRICS_FILE=$workspace/$file_name.metrics.txt AS=true VALIDATION_STRINGENCY=LENIENT

echo "$(date): running sambamba..."
## Sambamba
/shared/workspace/software/sambamba/0.4.7/bin/sambamba index -t 2 $workspace/$file_name.dedup.bam $workspace/$file_name.dedup.bam.bai


