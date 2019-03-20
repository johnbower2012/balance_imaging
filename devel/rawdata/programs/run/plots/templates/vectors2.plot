set multiplot layout 2,2;

set title 'eigenvector1';
vec=1
plot "eigenvectors0.dat" u vec title 'pipi' w lp, "eigenvectors1.dat" u vec title 'ppbar' w lp, "eigenvectors2.dat" u vec title 'pk' w lp, "eigenvectors3.dat" u vec title 'kk' w lp;

set title 'eigenvector2';
vec=2
plot "eigenvectors0.dat" u vec title 'pipi' w lp, "eigenvectors1.dat" u vec title 'ppbar' w lp, "eigenvectors2.dat" u vec title 'pk' w lp, "eigenvectors3.dat" u vec title 'kk' w lp;

title 'eigenvector3';
vec=3
plot "eigenvectors0.dat" u vec title 'pipi' w lp, "eigenvectors1.dat" u vec title 'ppbar' w lp, "eigenvectors2.dat" u vec title 'pk' w lp, "eigenvectors3.dat" u vec title 'kk' w lp;

set title 'eigenvector4';
vec=4
plot "eigenvectors0.dat" u vec title 'pipi' w lp, "eigenvectors1.dat" u vec title 'ppbar' w lp, "eigenvectors2.dat" u vec title 'pk' w lp, "eigenvectors3.dat" u vec title 'kk' w lp
