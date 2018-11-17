for((i=$1;i<$2;i++))
do 
fn=` printf "model_run%02d/" $i `
echo "Now running for directory ${fn} ..."
egrep -v "^('|$)" ${fn}mcmctrace.csv | tr -d ',' > ${fn}mcmctrace.dat
./extract.x 16000 4 20 ${fn}mcmctrace.dat ${fn}mcmcextract.dat
echo "Now ending for directory ${fn} ..."
echo ""
done
