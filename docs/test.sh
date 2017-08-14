#!/bin/bash

arr=(`cut -f 1- -d" " file.txt`)
for i in ${arr[@]}
do
    echo $i | sed -e "s/\n/ /g"
done
