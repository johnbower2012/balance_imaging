for((i=$1;i<$2;i++))
do
    fn=` printf "model_run%02d/" $i `
    echo "Running for directory ${fn} ..."
    cp -v calculations/moments/${fn}moments_plot.dat calculations/emulator/${fn}
    cp -v calculations/fitting/${fn}beta.dat calculations/emulator/${fn}
    echo "./calculations/emulator/emulator.x calculations/emulator/${fn}moments_plot.dat calculations/emulator/${fn}beta.dat calculations/emulator/${fn}hyperparameters.dat calculations/emulator/${fn}"
    ./calculations/emulator/emulator.x calculations/emulator/${fn}moments_plot.dat calculations/emulator/${fn}beta.dat calculations/emulator/${fn}hyperparameters.dat calculations/emulator/${fn}
    echo "train/test written to calculations/emulator/${fn}"
    echo "Finished running for directory ${fn} ..."
    echo ""
done