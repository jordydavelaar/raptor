
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

# Load prefix MultiNest files
prefix = sys.argv[1]

# Load data
MBHtot, M_UNITtot, Rhightot, itot, Chisquare_Radiotot = [[],[],[],[],[]]

MBH = np.loadtxt('%s_Chisquared_%d.txt'%(prefix,0),usecols=(0))
for i in range(len(MBH)):
	j = 0
	while True:
		try:
			data = np.loadtxt('%s_Chisquared_%d.txt'%(prefix,j))
		except:
			break
		MBH, M_UNIT, Rhigh, i2, Chisquare_Radio = data[i]
		MBHtot = np.append(MBHtot,MBH)
		M_UNITtot= np.append(M_UNITtot,M_UNIT)
		Rhightot = np.append(Rhightot,Rhigh)
		itot = np.append(itot,i2)
		Chisquare_Radiotot = np.append(Chisquare_Radiotot,Chisquare_Radio)
		j+=1
# Make figure
#print(MBHtot)
fig = plt.figure(num=None, figsize=(12,12))

# MBH plot
plt.subplot(321)

xlabelsize = '20'
ylabelsize = '20'
titlesize = '25'

plt.plot(MBHtot, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel(r'$MBH$ $M_{\odot}$', size = ylabelsize)
#plt.xlim
plt.ylim(1e7, 1e9)
#plt.xscale
plt.yscale('log')
plt.title('$MBH$ vs $Iteration$ $number$', size = titlesize)
#plt.legend(loc=0)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
#fig.savefig('M_UNITevolution.png')

# M_UNIT plot
plt.subplot(322)

xlabelsize = '20'
ylabelsize = '20'
titlesize = '25'

plt.plot(M_UNITtot, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$\mathcal{M}$ (g)', size = ylabelsize)
#plt.xlim
plt.ylim(1e22, 1e24)
#plt.xscale
plt.yscale('log')
plt.title('$\mathcal{M}$ vs $Iteration$ $number$', size = titlesize)
#plt.legend(loc=0)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
#fig.savefig('M_UNITevolution.png')



# Rhigh plot
#fig = plt.figure()
plt.subplot(323)

plt.plot(Rhightot, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$Rhigh$ ()', size = ylabelsize)
#plt.xlim
plt.ylim(1, 100)
#plt.xscale
#plt.yscale('lin')
plt.title('$Rhigh$ vs $Iteration$ $number$', size = titlesize)
#plt.legend(loc=0)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
#fig.savefig('Rhighevolution.png')



# i plot
#fig = plt.figure()
plt.subplot(324)

plt.plot(itot, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$i$ ($^\circ$)', size = ylabelsize)
#plt.xlim
plt.ylim(12, 45)
#plt.xscale
#plt.yscale('log')
plt.title('$i$ vs $Iteration$ $number$', size = titlesize)
#plt.legend(loc=0)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
#fig.savefig('ievolution.png')



# chi_squared plot
#fig = plt.figure()
plt.subplot(313)

Chisquare_Radio

plt.plot(Chisquare_Radiotot, color='b', label='Radio')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$\chi^2$ ()', size = ylabelsize)
#plt.xlim
#plt.ylim(12, 45)
#plt.xscale
plt.yscale('log')
plt.title('$\chi^2$ vs $Iteration$ $number$', size = titlesize)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
fig.savefig('%s_Convergence.png'%(prefix))







