#!/bin/bash

hello="goodbye.yo"

echo $hello | sed s/"\..*"//

echo ${hello//.*/}

echo `sed s/"\..*"// <<< $hello`

#hello="$hello/${$hello/\/.*}"

echo $hello
