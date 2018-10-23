./OBSERVABLES_EXP.sh calculations/data/exp_run00/ calculations/moments/exp_run00/
for((i=$1;i<$2;i++))
do
    fn=` printf "model_run%02d/" $i `
    ./OBSERVABLES_MODEL.sh calculations/data/${fn} calculations/moments/exp_run00/ calculations/moments/${fn}
done