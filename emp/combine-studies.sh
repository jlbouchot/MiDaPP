#!/usr/bin/env bash 

####
# Author: Gregory Ditzler
#         Drexel University
#         gregory.ditzler@gmail.com
#
# Merge map and biom files for the studies specified by the input 
# arguments. Ex.
#   >> ./combine-studies.sh /path/to/biom /path/to/map /new/file 1 2 3 4
# while look in /path/to/biom for the EMP biom files /path/to/maps
# for the map files. The studies indicated by 1, 2, 3 and 4 will be 
# merged. 

bpath=$1
shift 
mpath=$1
shift
npath=$1
shift

declare -A bfiles
declare -A mfiles
for arg in ${@}; do
  bfiles[$arg]=$bpath/study_${arg}_closed_reference_otu_table.biom
  mfiles[$arg]=$mpath/study_${arg}_mapping_file.txt
done
bioms=`echo ${bfiles[@]} | tr " " ","`
maps=`echo ${mfiles[@]} | tr " " ","`

# merge_mapping_files.py -m $maps -o $npath.txt
# merge_otu_tables.py -i $bioms -o $npath.biom
