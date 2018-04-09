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
markerstyles = [u'o',u's',u'^',r'$\square',r'$\circlearrowleft$',r'$\clubsuit$',r'$\checkmark$'
]
colors = ('r', 'g', 'b', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)

xx0=0.18
yy0=0.07
ww0=1.0-xx0-0.02
hh0=1.0-yy0-0.02

perfectdir='../../model_output/perfect/'
semiperfectdir='../../model_output/semiperfect/'
defaultdir='../../output_posterior_corrected/default/'
#modeldir='../../output_posterior_corrected/testrun/'
#modeldir='../../output_posterior_corrected/'+sys.argv[1]+'/'
#modeldir='../../output_posterior_corrected/onescale_1.0/'

for ipanel in range (0,4):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    filename=perfectdir+'I321_J2212.dat'
    perfectdata = np.loadtxt(filename,skiprows=6)
    filename=semiperfectdir+'I321_J2212.dat'
    semiperfectdata = np.loadtxt(filename,skiprows=6)
    filename=defaultdir+'I321_J2212.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pK.dat',skiprows=1)
    ymin=-0.05
    ymax=0.15
    ax.set_yticks(np.arange(-1,0.5,0.08), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.08), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(d) $K^-\!p$',family='serif', size=22,horizontalalignment='right')
  elif ipanel == 1:
    filename=perfectdir+'I2212_J2212.dat'
    perfectdata = np.loadtxt(filename,skiprows=6)
    filename=semiperfectdir+'I2212_J2212.dat'
    semiperfectdata = np.loadtxt(filename,skiprows=6)
    filename=defaultdir+'I2212_J2212.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_ppbar.dat',skiprows=1)
    ymin=-0.02
    ymax=0.45
    ax.set_yticks(np.arange(-0.2,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.6,0.05), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(c) $p\\bar{p}$',family='serif', size=22,horizontalalignment='right')
  elif ipanel == 2:
    filename=perfectdir+'I321_J321.dat'
    perfectdata = np.loadtxt(filename,skiprows=6)
    filename=semiperfectdir+'I321_J321.dat'
    semiperfectdata = np.loadtxt(filename,skiprows=6)
    filename=defaultdir+'I321_J321.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_KK.dat',skiprows=1)
    ymin=-0.02
    ymax=0.35
    ax.set_yticks(np.arange(0,0.6,0.1), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.1), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin),'(b) $K^+\!K^-$', family='serif', size=22,horizontalalignment='right')
  elif ipanel == 3:
    filename=perfectdir+'I211_J211.dat'
    perfectdata = np.loadtxt(filename,skiprows=6)
    filename=semiperfectdir+'I211_J211.dat'
    semiperfectdata = np.loadtxt(filename,skiprows=6)
    filename=defaultdir+'I211_J211.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pipi.dat',skiprows=1)
    ymin=-0.04
    ymax=1.0
    ax.set_yticks(np.arange(0,1.0,0.4), minor=False)
    ax.set_yticklabels(np.arange(0,1.0,0.4), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,1.0,0.1), minor=True)
    ax.text(1.7,ymin+0.8*(ymax-ymin), '(a) $\pi^+\!\pi^-$', family='serif', size=22,horizontalalignment='right')

  perfectdata = perfectdata.transpose()
  semiperfectdata = semiperfectdata.transpose()
  defaultdata = defaultdata.transpose()
  correcteddata = correcteddata.transpose()
  
  dely=defaultdata[0]
  B=defaultdata[1]
  dB=defaultdata[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=1,
  color=colors[0],markersize=8, marker=markerstyles[0],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor=colors[0])
  
  dely=perfectdata[0]
  B=perfectdata[1]
  dB=perfectdata[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=1,
  color=colors[1],markersize=8, marker=markerstyles[1],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor=colors[1])
  
  dely=semiperfectdata[0]
  B=semiperfectdata[1]
  dB=semiperfectdata[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=1,
  color=colors[2],markersize=8, marker=markerstyles[2],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor=colors[2])
  
  #dely=correcteddata[0]
  #B=correcteddata[1]
  #dB=correcteddata[2]
  #plt.errorbar(dely,B,yerr=dB,linestyle=linestyles[1],linewidth=1,
  #color=colors[0],markersize=8, marker=markerstyles[0],
  #markerfacecolor='none',
  #markeredgewidth=2,
  #markeredgecolor=colors[0])
  
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


plt.savefig('acceptance.pdf',format='pdf')
#plt.show()

os.system("open -a Preview acceptance.pdf")

quit()
