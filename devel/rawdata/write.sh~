samples=20
param=4
for((i=$1;i<$2;i++))
do
    mr=` printf "model_run_posterior%02d/" $i `
    echo "Running for directory ${mr} ..."
    for((j=0;j<$samples;j++))
    do
	dir=` printf "run%04d/" $j `
	cp -v ${mr}fixed_parameters.dat ${mr}${dir}
    done
    ./exec/write.x ${mr} mcmcextract.dat $samples $param
done