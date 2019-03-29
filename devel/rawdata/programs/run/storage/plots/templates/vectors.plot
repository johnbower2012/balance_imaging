set multiplot layout 2,2; 

set title 'eigenvector1'; 
plot "eigenvectors0.dat" u 1 title 'pipi' w lp, "eigenvectors1.dat" u 1 title 'ppbar' w lp, "eigenvectors2.dat" u 1 title 'pk' w lp, "eigenvectors3.dat" u 1 title 'kk' w lp; 

set title 'eigenvector2'; 
plot "eigenvectors0.dat" u 2 title 'pipi' w lp, "eigenvectors1.dat" u 2 title 'ppbar' w lp, "eigenvectors2.dat" u 2 title 'pk' w lp, "eigenvectors3.dat" u 2 title 'kk' w lp;

set title 'eigenvector3'; 
plot "eigenvectors0.dat" u 3 title 'pipi' w lp, "eigenvectors1.dat" u 3 title 'ppbar' w lp, "eigenvectors2.dat" u 3 title 'pk' w lp, "eigenvectors3.dat" u 3 title 'kk' w lp; 

set title 'eigenvector4'; 
plot "eigenvectors0.dat" u 4 title 'pipi' w lp, "eigenvectors1.dat" u 4 title 'ppbar' w lp, "eigenvectors2.dat" u 4 title 'pk' w lp, "eigenvectors3.dat" u 4 title 'kk' w lp