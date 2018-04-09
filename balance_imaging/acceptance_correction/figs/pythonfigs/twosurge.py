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

twosurgedir='../../output_posterior_corrected/default/'
onesurgedir='../../output_posterior_corrected/onescale_1.0/'
#twosurgedir='../../output_posterior_corrected/'+sys.argv[1]+'/'
#twosurgedir='../../output_posterior_corrected/onescale_1.0/'

for ipanel in range (0,4):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    filename=twosurgedir+'I321_J2212.dat'
    twosurgedata = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir+'I321_J2212.dat'
    onesurgedata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pK.dat',skiprows=1)
    ymin=-0.05
    ymax=0.13
    ax.set_yticks(np.arange(-1,0.5,0.08), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.08), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.annotate('$K^-p$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 1:
    filename=twosurgedir+'I2212_J2212.dat'
    twosurgedata = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir+'I2212_J2212.dat'
    onesurgedata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_ppbar.dat',skiprows=1)
    ymin=-0.02
    ymax=0.55
    ax.set_yticks(np.arange(-0.2,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.6,0.05), minor=True)
    ax.annotate(r'$p\bar{p}$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 2:
    filename=twosurgedir+'I321_J321.dat'
    twosurgedata = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir+'I321_J321.dat'
    onesurgedata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_KK.dat',skiprows=1)
    ymin=-0.02
    ymax=0.45
    ax.set_yticks(np.arange(0,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.annotate('$K^+K^-$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')
  elif ipanel == 3:
    filename=twosurgedir+'I211_J211.dat'
    twosurgedata = np.loadtxt(filename,skiprows=6)
    filename=onesurgedir+'I211_J211.dat'
    onesurgedata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pipi.dat',skiprows=1)
    ymin=-0.04
    ymax=1.0
    ax.set_yticks(np.arange(0,1.0,0.4), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.4), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=True)
    ax.annotate('$\pi^+\pi^-$',xy=(1.7,ymin+0.8*(ymax-ymin)), family='serif', size=24,horizontalalignment='right')

  twosurgedata = twosurgedata.transpose()
  onesurgedata = onesurgedata.transpose()
  correcteddata = correcteddata.transpose()

  dely=twosurgedata[0]
  B=twosurgedata[1]
  dB=twosurgedata[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=3,
  color=colors[0])

  dely=onesurgedata[0]
  B=onesurgedata[1]
  dB=onesurgedata[2]
  plt.plot(dely,B,linestyle=linestyles[2],linewidth=2,
  color=colors[1])
  
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


plt.savefig('twosurge.pdf',format='pdf')
#plt.show()

os.system("open -a Preview twosurge.pdf")

quit()
