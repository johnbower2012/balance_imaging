#!/bin/bash

folder="model_output"
filename[0]="I211_J211.dat"
filename[1]="I2212_J2212.dat"
filename[2]="I321_J2212.dat"
filename[3]="I321_J321.dat"
files=4
start=0
end=2

for((file=0;file<${files};file++))
do
    for((i=${start};i<${end};i++))
    do
	run=` printf "run%04d" $i `
	name=${folder}/${run}/${filename[${file}]}
	tempname=${folder}/${run}/temp
	printf "#" > ${tempname}
	cat ${name} >> ${tempname}
	cat ${tempname} > ${name}
	rm ${tempname}
    done
done