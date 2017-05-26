#!/bin/bash

file=$1

for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
grep "^SN" $file | grep "$stat" | awk -F $'\t' '{ print $2,"\t",$3}'
done


