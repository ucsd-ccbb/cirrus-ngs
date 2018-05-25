#!/bin/bash

pids=$(ps aux | grep Pipeline | grep -v stop_pipeline | tr -s ' ' | cut -f2 -d' ')

for pid in $pids
do
   kill $pid
done

qstat_jobs=$(qstat | grep ec2-user)

while read -r qstat_job
do
    qstat_job=$(echo $qstat_job | tr -s ' ' | cut -f1 -d' ')
    qdel $qstat_job
done <<< "$qstat_jobs"
