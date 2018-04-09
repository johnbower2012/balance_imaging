clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

read ./xyzt.dat tt xx y z px py

!lower panel
!___________________________________________
set xleadz 1.0
set yleadz 1.0
set xnumsz .4
set ynumsz .4
set xticl -.5
set xtics -.25
set yticl -.5
set ytics -.25
set xlabsz 1.0
set ylabsz 1.0
set lintyp 1.0
set linthk 6.0
set charsz 0.01
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
delxlabel=2.0
delylabel=2.5
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
set txthit 0.5
text `t (fm/c)'

! x-axis label
xtext=(xxuu+xxll)/2.0
ytext=yyll-delxlabel
set xloc xtext
set yloc ytext
set txtang 0
set txthit .5
text `x (fm)'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales -20 40.0 0 50

! These 8 par.s set the locations of tick marks
set xvmin -40.0
set xvmax 80.0
set nlxinc 6
set nsxinc 2

set yvmin 0
set yvmax 50
set nlyinc 5
set nsyinc 2

graph\axesonly


set pchar -12
color red
graph\noaxes xx tt

vector xa ta xb tb 2
xa(1)=2
ta(1)=0
xa(2)=xa(1)+0.2*50
ta(2)=50

!xb(1)=-2
!tb(1)=0
!xb(2)=xb(1)+0.9066*50
!tb(2)=50

color 1
set lintyp 5
set pchar 0
!graph\noaxes xa ta
!graph\noaxes xb tb

hardcopy s xyzt.eps

