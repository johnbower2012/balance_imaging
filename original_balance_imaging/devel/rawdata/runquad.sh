#! /bin/sh
make balance
ncores=4
nruns=1024
nsplit=`expr ${nruns} / ${ncores}`
echo nsplit=${nsplit}
for((i=0;i<${ncores};i++))
do
	firstrun=`expr ${i} \* ${nsplit}`
	lastrun=`expr ${firstrun} + ${nsplit} - 1`
	echo  firstrun=${firstrun}, lastrun=${lastrun}
	run.sh ${firstrun} ${lastrun} > logfiles/run${firstrun}_run${lastrun}.dat &
done