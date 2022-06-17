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
MBH, M_UNIT, Rhigh, i, Chisquare_Radio = np.loadtxt('%s_Chisquared.txt'%(prefix)).transpose()

# Make figure
fig = plt.figure(num=None, figsize=(12,12))

# MBH plot
plt.subplot(321)

xlabelsize = '20'
ylabelsize = '20'
titlesize = '25'

plt.plot(MBH, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$MBH$ (g)', size = ylabelsize)
#plt.xlim
plt.ylim(1e39, 1e41)
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

plt.plot(M_UNIT, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$\mathcal{M}$ (g)', size = ylabelsize)
#plt.xlim
plt.ylim(1e19, 1e24)
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

plt.plot(Rhigh, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$Rhigh$ ()', size = ylabelsize)
#plt.xlim
plt.ylim(1, 160)
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

plt.plot(i, '.')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$i$ ($^\circ$)', size = ylabelsize)
#plt.xlim
plt.ylim(0, 90)
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

plt.plot(Chisquare_Radio, color='b', label='Radio')
plt.xlabel('$Iteration$ $number$', size = xlabelsize)
plt.ylabel('$\chi^2$ ()', size = ylabelsize)
#plt.xlim
#plt.ylim(0, 90)
#plt.xscale
plt.yscale('log')
plt.title('$\chi^2$ vs $Iteration$ $number$', size = titlesize)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

fig.tight_layout()
fig.savefig('%s_Convergence.png'%(prefix))







