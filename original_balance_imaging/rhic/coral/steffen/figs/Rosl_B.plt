clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

read ../results/results3d_default.dat kt lambda_default Rout_default Routlab_default Rside_default Rlong_default ratio_default
read ../results/results3d_B0.dat kt lambda_B0 Rout_B0 Routlab_B0 Rside_B0 Rlong_B0 ratio_B0
read ../results/results3d_B4.dat kt lambda_B4 Rout_B4 Routlab_B4 Rside_B4 Rlong_B4 ratio_B4

read ../results/star.dat kt_star lambda Routlab_star Rside_star Rlong_star
ratio_star=Routlab_star/Rside_star

!__________________________________________
! LOWER PANEL
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
xxll=4.0
xxuu=16
yyll=3.5
yyuu=9
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=2.5
delylabel=2.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
set cursor -2
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `R<_>long<^>'

! x-axis label
xtext=(xxuu+xxll)/2.0
ytext=yyll-delxlabel
set xloc xtext
set yloc ytext
set txtang 0
set txthit .7
text `k<_>t<^> (MeV/c)'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 150 550 3 9.0

! These 8 par.s set the locations of tick marks
set xvmin 0.0
set xvmax 600
set nlxinc 6
set nsxinc 2

set yvmin 3.0
set yvmax 9.0
set nlyinc 3
set nsyinc 2

graph\axesonly

color green
set pchar 11
set lintyp 1
graph\noaxes kt Rlong_B0
color 1
set pchar 12
graph\noaxes kt Rlong_default
color blue
set pchar 10
graph\noaxes kt Rlong_B4

color red
set pchar -15
graph\noaxes kt_star Rlong_star


xtext=xxuu-1
set cursor -3
ytext=yyll+1.4
set xloc xtext
set yloc ytext
set txtang 0
set txthit 0.6
color green
!text `<eta>=2<eta,_>DSS<^> , <zeta><_>max<^>=2<eta,_>DSS<^>'

!-----------------------------------------------
!__________________________________________
! LOWER MIDDLE PANEL
!___________________________________________
set xnumsz 0.0
yyll=9
yyuu=14.5
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=2.8
delylabel=2.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
set cursor -2
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `R<_>side<^>'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 150 550 3.0001 7.0

! These 8 par.s set the locations of tick marks

set yvmin 3.0
set yvmax 9.0
set nlyinc 3
set nsyinc 2

graph\axesonly

color green
set pchar 11
set lintyp 1
graph\noaxes kt Rside_B0
color 1
set pchar 12
graph\noaxes kt Rside_default
color blue
set pchar 10
graph\noaxes kt Rside_B4

color red
set pchar -15
graph\noaxes kt_star Rside_star


!__________________________________________
! UPPER MIDDLE PANEL
!___________________________________________
set xnumsz 0.0
yyll=14.5
yyuu=20
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=2.8
delylabel=2.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
set cursor -2
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `R<_>out<^>'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 150 550 4 8.0

! These 8 par.s set the locations of tick marks

set yvmin 3.0
set yvmax 9.0
set nlyinc 3
set nsyinc 2

graph\axesonly

color green
set pchar 11
set lintyp 1
graph\noaxes kt Routlab_B0
color 1
set pchar 12
graph\noaxes kt Routlab_default
color blue
set pchar 10
graph\noaxes kt Routlab_B4

color red
set pchar -15
graph\noaxes kt_star Routlab_star


!-----------------------------------------------
!__________________________________________
! UPPER PANEL
!___________________________________________
set xnumsz 0.0
yyll=20
yyuu=24
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=2.8
delylabel=2.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
set cursor -2
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `R<_>out<^>/R<_>side<^>'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 150 550 0.8001 1.4

! These 8 par.s set the locations of tick marks

set yvmin 0.6
set yvmax 1.4
set nlyinc 2
set nsyinc 2

graph\axesonly

color green
set pchar 11
set lintyp 1
graph\noaxes kt ratio_B0
color 1
set pchar 12
graph\noaxes kt ratio_default
color blue
set pchar 10
graph\noaxes kt ratio_B4

color red
set pchar -15
graph\noaxes kt_star ratio_star


!-----------------------------------------------
hardcopy s Rosl_B.eps
