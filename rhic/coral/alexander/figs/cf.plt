clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

read ../results/cf_kt300.dat q cinv cout cside clong
read ../results/cf_gauss.dat q_gauss cinv_gauss cout_gauss cside_gauss clong_gauss

!lower panel
!___________________________________________
set xleadz 1.0
set yleadz 1.0
set xnumsz .6
set ynumsz .6
set xticl -.5
set xtics -.25
set yticl -.5
set ytics -.25
set xlabsz 1.0
set ylabsz 1.0
set lintyp 1.0
set linthk 6.0
set charsz 0.4
!set xlog 10
!set ylog 10
!!ylabel "P(n)"
!!xlabel "n"
!___________________________________________
!These are the dimensions of the borders
xxll=4.5
xxuu=15
yyll=3.0
yyuu=13
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=2.5
delylabel=3.0
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `C(q<_>inv<^>)-1'

! x-axis label
xtext=(xxuu+xxll)/2.0
ytext=yyll-delxlabel
set xloc xtext
set yloc ytext
set txtang 0
set txthit .7
text `q (MeV/c)'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 0 40.0 -0.1 0.7

! These 8 par.s set the locations of tick marks
set xvmin 0.0
set xvmax 100.0
set nlxinc 5
set nsxinc 4

set yvmin -1
set yvmax 1.0
set nlyinc 4
set nsyinc 5

graph\axesonly


set lintyp 1

set pchar 12
color red
graph\noaxes q cinv

color green
set pchar 10
graph\noaxes q_gauss cinv_gauss


hardcopy s cf.eps
