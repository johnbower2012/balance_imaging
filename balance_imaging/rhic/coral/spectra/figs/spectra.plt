clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

!___________________________________________
set xleadz 1.0
set yleadz 1.0
set xnumsz .5
set ynumsz .5
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
set ylog 10
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
delxlabel=2.6
delylabel=3.0
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
!label\yaxis `this is my label'
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.6
text `dN/(2<pi>p<_>t<^>dp<_>t<^>dy) (GeV/c)<^>-2<_>'

! x-axis label
xtext=(xxuu+xxll)/2.0
ytext=yyll-delxlabel
set xloc xtext
set yloc ytext
set txtang 0
set txthit .6
text `m<_>t<^>-m (GeV)'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales  0 2.0 -2 3

! These 8 par.s set the locations of tick marks
set xvmin 0.0
set xvmax 2.0
set nlxinc 4
set nsxinc 5

set yvmin 0.01
set yvmax 1000
set nlyinc 5
set nsyinc 10

graph\axesonly

read ../results/spectra_default.dat mt spi sk sp

scale=1
spi=spi*scale
!spi=spi/(1.0-exp(-(mt-0.07)/0.15))

scale=1
sk=sk*scale

set pchar 0
color red
graph\noaxes mt spi
color green
graph\noaxes mt sk
scale=0.7
sp=sp*scale
color cyan
!graph\noaxes mt sp

read ../results/spectra_noweak.dat mt spi sk sp

scale=0.7
sp=sp*scale

color blue
graph\noaxes mt sp

read phenixdata/phenix_piminus.dat pt_piminus mb_piminus mbe_piminus phenix_piminus phenix_piminus_error
mt_piminus=sqrt(pt_piminus*pt_piminus+0.13957*0.13957)-0.13957
set pchar -1
color red
graph\noaxes mt_piminus phenix_piminus phenix_piminus_error

read phenixdata/phenix_piplus.dat pt_piplus	 mb_piplus	 mbe_piplus	 phenix_piplus	 phenix_piplus_error
mt_piplus=sqrt(pt_piplus*pt_piplus+0.13957*0.13957)-0.13957
set pchar -12
color red
graph\noaxes mt_piplus	 phenix_piplus	phenix_piplus_error

read phenixdata/phenix_proton.dat pt_proton	 mb_proton	 mbe_proton	 phenix_proton	 phenix_proton_error
mt_proton=sqrt(pt_proton*pt_proton+0.9383*0.9383)-0.9383
set pchar -12
color blue
graph\noaxes mt_proton	 phenix_proton	phenix_proton_error

read phenixdata/phenix_pbar.dat pt_pbar	 mb_pbar	 mbe_pbar	 phenix_pbar	 phenix_pbar_error
mt_pbar=sqrt(pt_pbar*pt_pbar+0.9383*0.9383)-0.9383
set pchar -1
color blue
graph\noaxes mt_pbar	 phenix_pbar phenix_pbar_error

read phenixdata/phenix_kminus.dat pt_kminus mb_kminus mbe_kminus phenix_kminus phenix_kminus_error
mt_kminus=sqrt(pt_kminus*pt_kminus+0.494*0.494)-0.494
set pchar -1
color green
graph\noaxes mt_kminus phenix_kminus phenix_kminus_error

read phenixdata/phenix_kplus.dat pt_kplus	 mb_kplus	 mbe_kplus	 phenix_kplus	 phenix_kplus_error
mt_kplus=sqrt(pt_kplus*pt_kplus+0.494*0.494)-0.494
set pchar -12
color green
graph\noaxes mt_kplus	 phenix_kplus phenix_kplus_error

!read stardata/star_piminus.dat mt_piminus s_piminus e_piminus
!set pchar -12
!color red
!graph\noaxes mt_piminus s_piminus e_piminus

!read stardata/star_pbar.dat mt_pbar s_pbar e_pbar
!set pchar -12
!color blue
!graph\noaxes mt_pbar s_pbar e_pbar

!read stardata/star_proton.dat mt_proton s_proton e_proton
!set pchar -12
!color blue
!graph\noaxes mt_proton s_proton e_proton
!-----------------------------------------------


hardcopy s spectra.eps
