#!/bin/bash

## requirements
## install coreutils e.g. brew install coreutils (gcsplit, gsed, gsort)
## daniel restrepo-montoya

# split by ## marker - process each pairwise cluster-file - input: total.collinearity
for filename in ./*.collinearity
do
    gcsplit -z "${filename}" /##*/ {*}
done

# identify target proteins amont in clusters-files
for filename in ./*.list
do
    grep -f "${filename}" xx* > "${filename}".out
done

# parse by target - output target.collinearity
for filename in ./*.out
do
    gsed 's/:/\t/g' "${filename}" > "${filename}"_1
    gsort -t $'\t' -k 1,1 -V "${filename}"_1 > "${filename}"_2
    awk -v OFS='\t' '{print $1}' "${filename}"_2 > "${filename}"_3
    gsort -u -k1,1 -V "${filename}"_3 > "${filename}"_4
    { xargs cat < "${filename}"_4  ; } > "${filename}".collinearity
    rm xx* | rm *.list.out_*
done

