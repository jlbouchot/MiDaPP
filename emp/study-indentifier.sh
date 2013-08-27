#!/usr/bin/env bash 


# search the html file and create a mapping from the study numbers to the name of the 
# project title for each of the studies. you can download the emp.html file from 
# http://www.microbio.me/emp/
#
# just download the html webpage and name it emp.html and place it in the path of this
# bash script. 

grep '<a href="ftp://thebeast.colorado.edu/pub/QIIME_DB_Public_Studies/study' emp.html \
  | sed -e "s/.*Studies\/\(.*.*\)/\1/g" -e "s/\">/ /g" -e "s/<.*$//g" \
  | awk '{print $2,$1}' \
  > emp-descriptions.txt

