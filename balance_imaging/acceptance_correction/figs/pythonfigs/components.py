import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
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
  '_', '-', '--', ':','-.'
]
markerstyles = [
  u's',u'o',r'$\square',r'$\circlearrowleft$',r'$\clubsuit$',r'$\checkmark$'
]
colors = ('r', 'g', 'b', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)

xx0=0.18
yy0=0.07
ww0=1.0-xx0-0.02
hh0=1.0-yy0-0.02

datadir='../../output_posterior_corrected/default/'

for ipanel in range (0,4):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    filename=datadir+'I321_J2212.dat'
    bfdata = np.loadtxt(filename,skiprows=6)
    ymin=-0.1
    ymax=0.19
    ax.set_yticks(np.arange(-1,0.5,0.08), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.08), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin), '(d) $K^-\!p$',family='serif', size=22,ha='right')
  elif ipanel == 1:
    filename=datadir+'I2212_J2212.dat'
    bfdata = np.loadtxt(filename,skiprows=6)
    ymin=-0.15
    ymax=0.55
    ax.set_yticks(np.arange(-0.2,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.6,0.05), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(c) $p\\bar{p}$',family='serif', size=22,horizontalalignment='right')
  elif ipanel == 2:
    filename=datadir+'I321_J321.dat'
    bfdata = np.loadtxt(filename,skiprows=6)
    ymin=-0.05
    ymax=0.39
    ax.set_yticks(np.arange(0,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(b) $K^+\!K^-$',family='serif', size=22,horizontalalignment='right')
  elif ipanel == 3:
    filename=datadir+'I211_J211.dat'
    bfdata = np.loadtxt(filename,skiprows=6)
    ymin=-0.04
    ymax=1.0
    ax.set_yticks(np.arange(0,1.0,0.4), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.4), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(a) $\pi^+\!\pi^-$',family='serif', size=22,horizontalalignment='right')

  bfdata = bfdata.transpose()

  dely=bfdata[0]
  B=bfdata[1]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=3,color='k')
  
  dely=bfdata[0]
  B=bfdata[5]
  plt.plot(dely,B,linestyle=linestyles[2],linewidth=3,color=colors[0])
  
  dely=bfdata[0]
  B=bfdata[6]
  plt.plot(dely,B,linestyle=linestyles[3],linewidth=3,color=colors[1])
  
  dely=bfdata[0]
  B=bfdata[7]
  plt.plot(dely,B,linestyle=linestyles[4],linewidth=3,color=colors[2])
  
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


plt.savefig('components.pdf',format='pdf')
#plt.show()

os.system("open -a Preview components.pdf")

quit()
