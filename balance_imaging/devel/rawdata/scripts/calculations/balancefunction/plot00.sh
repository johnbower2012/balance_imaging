set macros;
set multiplot layout 4,2;
set key noautotitle; unset xlabel; unset ylabel;
set key font ",20";

XRANGE="0.1:1.8";
BF_OPTIONS="lc rgb 'blue' w lines"
DATA_OPTIONS="lc rgb 'red' lw 2 w points"
XTICS="set xtics (0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6); set xlabel '{/Symbol D}y' font ',20';"
NOXTICS="set xtics ('' 0.2,'' 0.4,'' 0.6,'' 0.8,'' 1.0,'' 1.2,'' 1.4,'' 1.6); unset xlabel;"
TMARGIN1="set tmargin at screen 0.90; set bmargin at screen 0.70;"
TMARGIN2="set tmargin at screen 0.70; set bmargin at screen 0.50;"
TMARGIN3="set tmargin at screen 0.50; set bmargin at screen 0.30;"
TMARGIN4="set tmargin at screen 0.30; set bmargin at screen 0.10;"
LMARGIN="set lmargin at screen 0.05; set rmargin at screen 0.5;"
RMARGIN="set lmargin at screen 0.50; set rmargin at screen 0.95;"

set title "Prior" font ",20";
set label "{/Symbol p}^+{/Symbol p}^-" enhanced font ",20" at 1.4,0.5;
YRANGE="-0.01:0.6"
YTICS="set ytics (0,0.1,0.2,0.3,0.4,0.5); unset ylabel;"
NOYTICS="set ytics ('' 0, '' 0.1, '' 0.2, '' 0.3, '' 0.4, '' 0.5); unset ylabel;"
@NOXTICS; @YTICS;
@TMARGIN1; @LMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/prior/pipi_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:1002 @DATA_OPTIONS;
set title "Posterior" font ",20";
@NOXTICS; @NOYTICS;
@TMARGIN1; @RMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/posterior/pipi_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:22 @DATA_OPTIONS;
unset label;
unset title;

set label "p ~p{.4-}" enhanced font ",20" at 1.4,0.3;
set label "B({/Symbol D}y)" enhanced font ",20" at -0.05,-0.05 rotate left;
YRANGE="-0.01:0.4"
YTICS="set ytics (0,0.05,0.1,0.15,0.2,0.25,0.3,0.35); unset ylabel;"
NOYTICS="set ytics ('' 0,'' 0.05, '' 0.1, '' 0.15, '' 0.2, '' 0.25, '' 0.3, ''0.35); unset ylabel;"
@NOXTICS; @YTICS;
@TMARGIN2; @LMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/prior/ppbar_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:1002 @DATA_OPTIONS;
@NOXTICS; @NOYTICS; 
@TMARGIN2; @RMARGIN;
unset label;
set label "p ~p{.4-}" enhanced font ",20" at 1.4,0.3;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/posterior/ppbar_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:22 @DATA_OPTIONS;
unset label;

set label  "K^+K^-" enhanced font ",20" at 1.4,0.2;
YRANGE="-0.01:0.25"
YTICS="set ytics (0,0.05,0.1,0.15,0.2); unset ylabel;"
NOYTICS="set ytics ('' 0,'' 0.05, '' 0.1, '' 0.15, '' 0.2); unset ylabel;"
@TMARGIN3; @LMARGIN;
@NOXTICS; @YTICS;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/prior/kk_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:1002 @DATA_OPTIONS;
@NOXTICS; @NOYTICS;
@TMARGIN3; @RMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/posterior/kk_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:22 @DATA_OPTIONS;
unset label;

set label  "K^-p" enhanced font ",20" at 1.4,0.06;
YRANGE="-0.03:0.08"
YTICS="set ytics (-0.02,0.0,0.02,0.04,0.06); unset ylabel;"
NOYTICS="set ytics ('' -0.02,'' 0.0,'' 0.02,'' 0.04,'' 0.06); unset ylabel;"
@XTICS; @YTICS;
@TMARGIN4; @LMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/prior/pk_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:1002 @DATA_OPTIONS;
@XTICS; @NOYTICS;
@TMARGIN4; @RMARGIN;
plot [@XRANGE][@YRANGE] for[i=2:21] "model_run00/posterior/pk_exp.dat" every ::1::17 u 1:i @BF_OPTIONS, "" every ::1::17 u 1:22 @DATA_OPTIONS;
unset label;