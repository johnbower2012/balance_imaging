for((i=0;i<5;i++))
do
    fn=` printf "model_run%02d/" $i`
    echo Now copying for ${fn} ...
    cp -v calculations/moments/${fn}moments_plot.dat calculations/mcmc/${fn}
    cp -v calculations/fitting/${fn}beta.dat calculations/mcmc/${fn}
    cp -v calculations/mcmc/hyperparameters.dat calculations/mcmc/${fn}
    cp -v calculations/moments/${fn}moments_exp_z.dat calculations/mcmc/${fn}
done