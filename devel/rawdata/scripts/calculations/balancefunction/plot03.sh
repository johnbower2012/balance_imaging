set multiplot layout 4,2 title "Model Balance Function BF=0.25"; 
set key noautotitle; set grid; set xlabel "Relative Rapidity"; 
set ylabel "PiPi Balance Function"; 
set title "Prior";plot [0:1.8][-0.2:0.9] for[i=2:21] "model_run03/prior/pipi_exp.dat" u 1:i w lines, "" u 1:1002 lc rgb "black" lw 3 w lines title "experimental BF"; 
set title "Posterior"; plot [0:1.8][-0.2:0.9] for[i=2:21] "model_run03/posterior/pipi_exp.dat" u 1:i w lines, "" u 1:22 lc rgb "black" lw 3 w lines title "experimental BF"; 
set ylabel "PPbar Balance Function"; 
set title "Prior"; plot [0:1.8][-0.01:0.6] for[i=2:21] "model_run03/prior/ppbar_exp.dat" u 1:i w lines, "" u 1:1002 lc rgb "black" lw 3 w lines title "experimental BF"; 
set title "Posterior"; plot [0:1.8][-0.01:0.6] for[i=2:21] "model_run03/posterior/ppbar_exp.dat" u 1:i w lines, "" u 1:22 lc rgb "black" lw 3 w lines title "experimental BF"; 
set ylabel "KK Balance Function"; 
set title "Prior"; plot [0:1.8][-0.15:0.40] for[i=2:21] "model_run03/prior/kk_exp.dat" u 1:i w lines, "" u 1:1002 lc rgb "black" lw 3 w lines title "experimental BF"; 
set title "Posterior"; plot [0:1.8][-0.15:0.40] for[i=2:21] "model_run03/posterior/kk_exp.dat" u 1:i w lines, "" u 1:22 lc rgb "black" lw 3 w lines title "experimental BF"; 
set ylabel "pK Balance Function"; 
set title "Prior"; plot [0:1.8][-0.05:0.1] for[i=2:21] "model_run03/prior/pk_exp.dat" u 1:i w lines, "" u 1:1002 lc rgb "black" lw 3 w lines title "experimental BF"; 
set title "Posterior"; plot [0:1.8][-0.05:0.1] for[i=2:21] "model_run03/posterior/pk_exp.dat" u 1:i w lines, "" u 1:22 lc rgb "black" lw 3 w lines title "experimental BF";