clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

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
xxuu=17
yyll=5.0
yyuu=15
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=3
delylabel=3.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
label\yaxis `<delta>'
!ytext=(yyll+yyuu)/2.0
!xtext=xxll-delylabel
!set xloc xtext
!set yloc ytext
!set txtang 90
!set txthit 0.7
!text `dN/(p<_>t<^>dp<_>t<^>dy) arb. units'

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
scales 0 350.0 -20 20

! These 8 par.s set the locations of tick marks
set xvmin 0.0
set xvmax 400.0
set nlxinc 4
set nsxinc 4

set yvmin -20
set yvmax 20
set nlyinc 4
set nsyinc 5

graph\axesonly

read pk_phaseshifts.dat q d_s11 dddq_s11 d_p11 dddq_p11 d_p13 dddq_p13
d_s11=d_s11*(180.0/pi)
d_p11=d_p11*(180.0/pi)
d_p13=d_p13*(180.0/pi)
dddq_s11=dddq_s11*100*(180.0/pi)
dddq_p11=dddq_p11*100*(180.0/pi)
dddq_p13=dddq_p13*100*(180.0/pi)


set lintyp 1
set pchar 0
color red
graph\noaxes q d_s11
color green
graph\noaxes q d_p11
color blue
graph\noaxes q d_p13

set lintyp 3
set charsz 0.1
set pchar 0
color red
graph\noaxes q dddq_s11
color green
graph\noaxes q dddq_p11
color blue
graph\noaxes q dddq_p13

!-----------------------------------------------

hardcopy s pk_phaseshifts.eps
