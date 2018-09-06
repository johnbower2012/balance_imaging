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

Npanels=3

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



optionBdir='../../output_posterior_corrected/sigmaB_0.1/'
optionAdir='../../output_posterior_corrected/sigmaB_0.7/'
defaultdir='../../output_posterior_corrected/default/'
#optionBdir='../../output_posterior_corrected/'+sys.argv[1]+'/'
#optionBdir='../../output_posterior_corrected/onescale_1.0/'


for ipanel in range (0,Npanels):
  ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
  if ipanel == 0:
    optionBdir='../../output_posterior_corrected/sigmaB_0.1/'
    optionAdir='../../output_posterior_corrected/sigmaB_0.7/'
    filename=defaultdir+'I321_J2212.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    filename=optionBdir+'I321_J2212.dat'
    optionBdata = np.loadtxt(filename,skiprows=6)
    filename=optionAdir+'I321_J2212.dat'
    optionAdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_pK.dat',skiprows=1)
    ymin=-0.08
    ymax=0.16
    ax.set_yticks(np.arange(-1,0.5,0.08), minor=False)
    ax.set_yticklabels(np.arange(-1,0.5,0.08), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-1,0.5,0.02), minor=True)
    ax.text(1.4,ymin+0.8*(ymax-ymin),'(c) $K^-\!p$',family='serif', size=22,horizontalalignment='right')
    labelA='$\sigma_B/\sigma_A=0.7$'
    labelB='$\sigma_B/\sigma_A=0.1$'
    defaultlabel='$\sigma_B/\sigma_A=0.4$'
  elif ipanel == 1:
    optionBdir='../../output_posterior_corrected/QoverS_0.1/'
    optionAdir='../../output_posterior_corrected/QoverS_0.25/'
    filename=defaultdir+'I2212_J2212.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    filename=optionBdir+'I2212_J2212.dat'
    optionBdata = np.loadtxt(filename,skiprows=6)
    filename=optionAdir+'I2212_J2212.dat'
    optionAdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_ppbar.dat',skiprows=1)
    ymin=-0.02
    ymax=0.6
    ax.set_yticks(np.arange(-0.2,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(-0.2,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(-0.2,0.6,0.05), minor=True)
    ax.text(1.4,ymin+0.1*(ymax-ymin),'(b) $p\\bar{p}$',family='serif', size=22,horizontalalignment='right')
    labelA='Tr$\chi/s=0.25$'
    defaultlabel='Tr$\chi/s=0.18$'
    labelB='Tr$\chi/s=0.1$'
  elif ipanel == 2:
    optionBdir='../../output_posterior_corrected/SoverU_0.5/'
    optionAdir='../../output_posterior_corrected/SoverU_1.4/'
    filename=defaultdir+'I321_J321.dat'
    defaultdata = np.loadtxt(filename,skiprows=6)
    filename=optionBdir+'I321_J321.dat'
    optionBdata = np.loadtxt(filename,skiprows=6)
    filename=optionAdir+'I321_J321.dat'
    optionAdata = np.loadtxt(filename,skiprows=6)
    correcteddata = np.loadtxt('../../exp_data_corrected/star_KK.dat',skiprows=1)
    ymin=-0.02
    ymax=0.45
    ax.set_yticks(np.arange(0,0.6,0.2), minor=False)
    ax.set_yticklabels(np.arange(0,0.6,0.2), minor=False, family='serif', size=14)
    ax.set_yticks(np.arange(0,0.6,0.05), minor=True)
    ax.text(1.4,ymin+0.8*(ymax-ymin),'(a) $K^+\!K^-$',family='serif', size=22,horizontalalignment='right')
    labelB='$\chi_{ss}/\chi_{uu}=0.5$'
    labelA='$\chi_{ss}/\chi_{uu}=1.4$'
    defaultlabel='$\chi_{ss}/\chi_{uu}=0.93$'

  defaultdata = defaultdata.transpose()
  optionBdata = optionBdata.transpose()
  optionAdata = optionAdata.transpose()
  correcteddata = correcteddata.transpose()

  dely=optionBdata[0]
  B=optionBdata[1]
  dB=optionBdata[2]
  plt.plot(dely,B,linestyle=linestyles[2],linewidth=3,
  color=colors[1],label=labelB)

  dely=defaultdata[0]
  B=defaultdata[1]
  dB=defaultdata[2]
  plt.plot(dely,B,linestyle=linestyles[1],linewidth=3,
  color=colors[0],label=defaultlabel)
  
  dely=optionAdata[0]
  B=optionAdata[1]
  dB=optionAdata[2]
  plt.plot(dely,B,linestyle=linestyles[3],linewidth=2,
  color=colors[2],label=labelA)
  
  dely=correcteddata[0]
  B=correcteddata[1]
  dB=correcteddata[2]
  plt.errorbar(dely,B,yerr=dB,linestyle=linestyles[0],linewidth=2,
  color='k',markersize=8, marker=markerstyles[0],
  markerfacecolor='none',
  markeredgewidth=2,
  markeredgecolor='k')
  
  if ipanel == 0:
    legend(loc=3,fontsize=16)
  if ipanel == 1:
    legend(loc=(0.25,0.03),fontsize=16)
  if ipanel == 2:
    legend(loc=3,fontsize=16)
  
  ax.set_xticks(np.arange(0,2,0.5), minor=False)
  ax.set_xticks(np.arange(0,2,0.1), minor=True)
  if ipanel == 0:
    ax.set_xticklabels(np.arange(0,2.0,0.5), minor=False, size=14)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    plt.xlabel('$\Delta y$',fontsize=24)
    plt.ylabel('$B(\Delta y)$',fontsize=24,y=1.5)
    #plt.title('Balance Functions',fontsize=12, color='gray')
  else:
    ax.set_xticklabels([], minor=False)
  plt.xlim(0.0,1.49)

  plt.ylim(ymin,ymax)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(sformatter)


plt.savefig('parsens.pdf',format='pdf')
#plt.show()

os.system("open -a Preview parsens.pdf")

quit()
