#!/bin/bash
#this is to be used to plot

list=`find run*/$1`

for file in $list
do
    (egrep -v "^(#|$)" $file | tr -s ' ' | cut -d ' ' -f 1,2) >> $1
done
