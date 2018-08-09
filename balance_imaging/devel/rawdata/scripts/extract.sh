SOURCE_FOLDER='calculations/mcmc/'

for((i=$1;i<$2;i++))
do
    fn=` printf "model_run%02d/" $i`
    dn=` printf "model_run_posterior%02d/" $i `
    echo Now running for directory ${fn}...
    tail -16000 ${SOURCE_FOLDER}${fn}mcmctrace.csv | tr -d ',' > ${SOURCE_FOLDER}${fn}mcmctrace.dat
    ./${SOURCE_FOLDER}extract.x 16000 4 20 ${SOURCE_FOLDER}${fn}mcmctrace.dat ${SOURCE_FOLDER}${fn}mcmcextract.dat
    cp -v ${SOURCE_FOLDER}${fn}mcmcextract.dat ../${dn}
    echo "Ending directory ${fn}..."
    echo ""
done