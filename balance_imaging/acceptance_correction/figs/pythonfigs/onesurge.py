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
  '_', '-', '--', ':'
]
markerstyles = [
  u's',
  u'o',
  r'$\square',
  r'$\circlearrowleft$',
  r'$\clubsuit$',
  r'$\checkmark$'
]
colors = ('r', 'g', 'b', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)

xx0=0.18
yy0=0.07
ww0=1.0-xx0-0.02
hh0=1.0-yy0-0.02

onesurgedir_0='../../output_posterior_corrected/onescale_0.0/'
onesurgedir_5='../../output_posterior_corrected/onescale_0.5/'
onesurgedir_10='../../output_posterior_corrected/onescale_1.0/'


for ipanel in range (0,4):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    bfile='I321_J2212.dat'
    filename=onesurgedir_0+bfile
    data_0 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_5+bfile
    data_5 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_10+bfile
    data_10 = np.loadtxt(filename,skiprows=6)    
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pK.dat',skiprows=1)
    ymin=-0.05
    ymax=0.13
    ax.set_yticks(np.arange(-1,0.5,0.08), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.08), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.annotate('$K^-p$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 1:
    bfile='I2212_J2212.dat'
    filename=onesurgedir_0+bfile
    data_0 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_5+bfile
    data_5 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_10+bfile
    data_10 = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_ppbar.dat',skiprows=1)
    ymin=-0.02
    ymax=0.75
    ax.set_yticks(np.arange(-0.2,0.8,0.2), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.8,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.8,0.1), minor=True)
    ax.annotate(r'$p\bar{p}$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 2:
    bfile='I321_J321.dat'
    filename=onesurgedir_0+bfile
    data_0 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_5+bfile
    data_5 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_10+bfile
    data_10 = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_KK.dat',skiprows=1)
    ymin=-0.02
    ymax=0.45
    ax.set_yticks(np.arange(0,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.annotate('$K^+K^-$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 3:
    bfile='I211_J211.dat'
    filename=onesurgedir_0+bfile
    data_0 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_5+bfile
    data_5 = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir_10+bfile
    data_10 = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pipi.dat',skiprows=1)
    ymin=-0.04
    ymax=1.0
    ax.set_yticks(np.arange(0,1.0,0.4), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.4), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=True)
    ax.annotate('$\pi^+\pi^-$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  data_0 = data_0.transpose()
  data_5 = data_5.transpose()
  data_10 = data_10.transpose()
  correcteddata = correcteddata.transpose()

  dely=data_10[0]
  B=data_10[1]
  dB=data_10[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=3,
  color=colors[0])
  
  dely=data_5[0]
  B=data_5[1]
  dB=data_5[2]
  plt.plot(dely,B,linestyle=linestyles[2],linewidth=3,
  color=colors[1])
  
  dely=data_0[0]
  B=data_0[1]
  dB=data_0[2]
  plt.plot(dely,B,linestyle=linestyles[3],linewidth=3,
  color=colors[2])
  
  dely=correcteddata[0]
  B=correcteddata[1]
  dB=correcteddata[2]
  plt.errorbar(dely,B,yerr=dB,linestyle=linestyles[0],linewidth=2,
  color='k',markersize=8, marker=markerstyles[0],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor='k')
  
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


plt.savefig('onesurge.pdf',format='pdf')
#plt.show()

os.system("open -a Preview onesurge.pdf")

quit()
