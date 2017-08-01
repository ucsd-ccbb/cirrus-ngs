#!/bin/sh

mkdir $inputFolder"/output"

cat $inputFolder"/alignment"*".e"* > $outputFile

rm $inputFolder"/alignment"*"."*

echo "done."
