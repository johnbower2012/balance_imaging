import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import numpy as np
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

Npanels=4

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)
plt.rc('text', usetex=False)
linestyles = [
  '_', '-', '--', ':'
]
markerstyles = [
  u's', u'o', r'$\square', r'$\circlearrowleft$', r'$\clubsuit$', r'$\checkmark$'
]
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)

xx0=0.18
yy0=0.07
ww0=1.0-xx0-0.02
hh0=1.0-yy0-0.02

for ipanel in range (0,4):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    expdata = np.loadtxt('../../../exp_data/star_pK.dat',skiprows=1)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pK.dat',skiprows=1)
    ymin=-0.05
    ymax=0.1
    ax.set_yticks(np.arange(-1,0.5,0.04), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.04), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.text(1.75,0.07,'(d) $K^-\!p$',family='serif', size=20, ha='right')
  elif ipanel == 1:
    expdata = np.loadtxt('../../../exp_data/star_ppbar.dat',skiprows=1)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_ppbar.dat',skiprows=1)
    ymin=-0.02
    ymax=0.39
    ax.set_yticks(np.arange(-0.2,0.6,0.1), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.6,0.1), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.6,0.05), minor=True)
    ax.text(1.75,0.3,'(c) $p\\bar{p}$',family='serif', size=20, ha='right')
  elif ipanel == 2:
    expdata = np.loadtxt('../../../exp_data/star_KK.dat',skiprows=1)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_KK.dat',skiprows=1)
    ymin=-0.02
    ymax=0.35
    ax.set_yticks(np.arange(0,0.6,0.1), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.1), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.text(1.75,0.285,'(b) $K^+\!K^-$', family='serif', size=20, ha='right')
    ax.text(0.2,0.1,"STAR PRELIMINARY",family='serif',color=colors[1],size=18)
    
  elif ipanel == 3:
    expdata = np.loadtxt('../../../exp_data/star_pipi.dat',skiprows=1)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pipi.dat',skiprows=1)
    ymin=-0.04
    ymax=0.9
    ax.set_yticks(np.arange(0,1.0,0.2), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=True)
    ax.text(1.75,0.75,'(a) $\pi^+\pi^-$',family='serif', size=20, ha='right')

  expdata = expdata.transpose()
  correcteddata = correcteddata.transpose()
  dely=expdata[0]
  B=expdata[1]
  dB=expdata[2]

  plt.errorbar(dely,B,yerr=dB,linestyle=linestyles[1],linewidth=1,
  color=colors[1],markersize=8, marker=markerstyles[1],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor=colors[1])
  
  dely=correcteddata[0]
  B=correcteddata[1]
  dB=correcteddata[2]
  
  plt.errorbar(dely,B,yerr=dB,linestyle=linestyles[1],linewidth=1,
  color=colors[2],markersize=8, marker=markerstyles[0],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor=colors[2])
  
  ax.set_xticks(np.arange(0,2.0,0.5), minor=False)
  ax.set_xticks(np.arange(0,2.0,0.1), minor=True)
  if ipanel == 0:
    ax.set_xticklabels(np.arange(0,2.0,0.5), minor=False, size=14)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    plt.xlabel('$\Delta y$',fontsize=24)
    plt.ylabel('$B(\Delta y)$',fontsize=24,y=2.0)
    #plt.title('Balance Functions',fontsize=12, color='gray')
  else:
    ax.set_xticklabels([], minor=False)
  plt.xlim(0.0,1.8)

  plt.ylim(ymin,ymax)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(sformatter)


  plt.savefig('Bcorrected_STAR.pdf',format='pdf')
#plt.show()

os.system('open -a Preview Bcorrected_STAR.pdf')

quit()
