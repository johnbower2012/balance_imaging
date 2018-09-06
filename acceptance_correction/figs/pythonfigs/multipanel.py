import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
linestyles = [
  '_', '-', '--', ':'
]
markerstyles = [
  r'$\lambda$',
  r'$\bowtie$',
  r'$\circlearrowleft$',
  r'$\clubsuit$',
  r'$\checkmark$'
]
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.5,0.8,0.4])

#x = np.arange(0,40.1,1)
#y = x**2

mydata = np.loadtxt('junk.dat',skiprows=1)
mydata = mydata.transpose()
x=80.0*mydata[0]
y=15000.0*mydata[1]
z=10000*mydata[2]

plt.plot(x,y,linestyle=linestyles[1],linewidth=2,color=colors[1],markersize=8, marker=markerstyles[1], markerfacecolor=colors[1], markeredgecolor=colors[1])
#plt.semilogy(x,y)
ax.set_xticks(np.arange(0,41,10), minor=False)
ax.set_xticklabels([], minor=False)
ax.set_xticks(np.arange(0,41,5), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

ax.set_yticks(np.arange(0,200000,40000), minor=False)
ax.set_yticklabels(np.arange(0,200000,40000), minor=False, family='serif')
ax.set_yticks(np.arange(0,200000,10000), minor=True)
plt.ylim(0.0,160000.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

#plt.xlabel('time (s)')
plt.ylabel('voltage=$RI$ (mV)',fontsize=16)
plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
fontsize=12, color='gray')
# Make room for the ridiculously large title.
plt.subplots_adjust(top=0.85)

# LOWER PANEL

ax = fig.add_axes([0.15,0.1,0.8,0.4])

plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color=colors[3],markersize=8, marker=markerstyles[3], markerfacecolor=colors[3], markeredgecolor=colors[3])
#plt.semilogy(x,y)
ax.set_xticks(np.arange(0,41,10), minor=False)
ax.set_xticklabels(np.arange(0,41,10), minor=False)
ax.set_xticks(np.arange(0,41,5), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

ax.set_yticks(np.arange(0,200000,40000), minor=False)
ax.set_yticklabels(np.arange(0,200000,40000), minor=False, family='serif')
ax.set_yticks(np.arange(0,200000,10000), minor=True)
plt.ylim(0.0,159500.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('time (s)')
plt.ylabel('voltage=$RI$ (mV)',fontsize=16)
plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
fontsize=12, color='gray')
# Make room for the ridiculously large title.
plt.subplots_adjust(top=0.85)


plt.savefig('multipanel.pdf',format='pdf')
#plt.show()

quit()
