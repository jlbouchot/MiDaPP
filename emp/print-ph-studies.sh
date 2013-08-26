#!/usr/bin/env bash

cpath=`pwd`
for f in `find maps/study*.txt`; do 
  g=`python $cpath/../map_analyzer.py --map $f --print_meta | sed -e "s/ //g" |grep "PH$"`
  if [ "$g" == "PH" ]; then 
    echo "$f" 
  fi 
done

