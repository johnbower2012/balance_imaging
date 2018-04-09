#! /bin/sh
for((i=$1;i<=$2;i++))
do
	run_name=`printf "run%04d" ${i}`;
	echo  XXXXXXXX  run_name=${run_name}  XXXXXXXX;
	balance ${run_name};
done