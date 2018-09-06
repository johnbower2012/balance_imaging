for((i=$1;i<$2;i++))
do
    fn=` printf "model_run%02d/" $i`
    echo Now running for directory ${fn}...
    echo "./calculations/fitting/fitting.x calculations/moments/${fn}moments_plot.dat calculations/fitting/${fn}"
    ./calculations/fitting/fitting.x calculations/moments/${fn}moments_plot.dat calculations/fitting/${fn}
    cp -v calculations/fitting/${fn}beta.dat calculations/mcmc/${fn}
    echo "Ending directory ${fn}..."
    echo ""
done