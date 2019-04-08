set terminal pngcairo size 1400, 1200
set output 'stacking.png'
set lmargin at screen 0.08
set rmargin at screen 0.95
TOP=0.95
DY=0.22
start=2; finish=201; inc=40;
set multiplot layout 4,1;
set xrange[0.15:1.5]
set xtics font ",24"
set ytics font ",24"
set offset 0,0,graph 0.05, graph 0.05
unset key;

set xlabel '{/Symbol D}Y, Rapidity' offset 0,-2 enhanced font ",26"
set ylabel 'KK' offset -3 font ",32"
set ytics ("-0.1" -0.1,"0" 0,"0.1" 0.1,"0.2" 0.2,"0.3" 0.3)
set tmargin at screen TOP-3*DY
set bmargin at screen TOP-4*DY
oname="I3.dat"
rname="recon3.dat"
plot for[i=start:finish:inc] oname u 1:i lc rgb 'black' lw 3 w l, for[i=start:finish:inc] rname u 1:i lw 3 lc rgb 'blue' w l

unset xlabel
set xtics format ''
set ytics ("-0.1" -0.1,"" -0.05, "0" 0,"" 0.05,"0.1" 0.1)
set ylabel 'pK' offset -3 font ",32"
set tmargin at screen TOP-2*DY
set bmargin at screen TOP-3*DY
oname="I2.dat"
rname="recon2.dat"
plot for[i=start:finish:inc] oname u 1:i lc rgb 'black' lw 3 w l, for[i=start:finish:inc] rname u 1:i lw 3 lc rgb 'blue' w l

set ylabel 'p~p{.4-}' enhanced offset -4 font ",32"
set ytics ("" -0.1,"0" 0,"0.1" 0.1,"0.2" 0.2,"0.3" 0.3,"0.4" 0.4)
set tmargin at screen TOP-1*DY
set bmargin at screen TOP-2*DY
oname="I1.dat"
rname="recon1.dat"
plot for[i=start:finish:inc] oname u 1:i lw 3 lc rgb 'black' w l, for[i=start:finish:inc] rname u 1:i lw 3 lc rgb 'blue' w l

set title 'BF Reconstruction keeping 10 Princ Comp' font ",32";
set ylabel '{/Symbol p}{/Symbol p}' enhanced offset -4 font ",32"
set ytics ("-0.1" -0.1,"" 0,"0.1" 0.1,"" 0.2,"0.3" 0.3,"" 0.4, "0.5" 0.5)
set tmargin at screen TOP-0*DY
set bmargin at screen TOP-1*DY
oname="I0.dat"
rname="recon0.dat"
set label "Model" textcolor rgb "black" font ",32" at 1.2,0.36
set label "Reconstruction" textcolor rgb "blue" font ",32" at 1.2,0.30
plot for[i=start:finish:inc] oname u 1:i lw 3 lc rgb 'black' w l, for[i=start:finish:inc] rname u 1:i lw 3 lc rgb 'blue' w l

unset multiplot; set output
