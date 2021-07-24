#!/bin/bash
#For every file in results Directory
for d in */ ; do

#Enter into Directory
cd $d

name=$(echo $d | sed 's:/[^/]*$::')

echo $name

awk 'NR >13 { print }' < consurf.grades > $name.clean.grade

cd ..

done
