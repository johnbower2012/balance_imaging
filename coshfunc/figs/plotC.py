import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.18,0.12,0.8,0.8])

mydata = np.loadtxt('crap.dat',skiprows=0,unpack=True)
eta=mydata[0]
C0=mydata[1]
C1=mydata[2]
C2=mydata[3]
C3=mydata[4]
C4=mydata[5]
C5=mydata[6]
C6=mydata[7]
C7=mydata[8]
C8=mydata[9]

#x = np.arange(0,40.1,1)
#y =100.0* x**2
#z=0.5*y
#Use linestyle='None' for no line...
plt.plot(eta,C0,linestyle='-',color='r')
plt.plot(eta,C1,linestyle='-',color='b')
plt.plot(eta,C2,linestyle='-',color='g')
plt.plot(eta,C3,linestyle='-',color='y')
plt.plot(eta,C4,linestyle='-',color='cyan')
#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,10,1), minor=False)
ax.set_xticklabels(np.arange(0,10,1), minor=False, family='serif')
ax.set_xticks(np.arange(0,10,0.5), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,6.0)

ax.set_yticks(np.arange(-1,2,0.25), minor=False)
ax.set_yticklabels(np.arange(-1,2.0,0.25), minor=False, family='serif')
ax.set_yticks(np.arange(-1,2,0.05), minor=True)
plt.ylim(-1.05,1.05)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\eta$', fontsize=18, weight='normal')
plt.ylabel('$C(\eta)$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('plotC.pdf',format='pdf')
os.system('open -a Preview plotC.pdf')
#plt.show()
quit()

