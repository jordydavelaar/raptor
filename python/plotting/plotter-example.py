#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt
import h5py
import rapplot

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

ind = int(sys.argv[1])

#need to know the data_ids
data_id = rapplot.read_data_id("output",ind)

#Computing relevant constants
M=5e6 * rapplot.MSUN
d=8 * rapplot.KPC

rg = (rapplot.G*M/rapplot.SPEED_OF_LIGHT**2.)

mas = (rg/d)* rapplot.MAS_IN_DEG

Tunit =rg/rapplot.SPEED_OF_LIGHT

halfrange=20 #in rg

#read data
min,max,image = rapplot.read_data("output",ind,data_id)

#plot data
plt.figure(figsize=(6,5),dpi=1000,facecolor='w')
fig, axs = plt.subplots(2,2,figsize=(12,10))

stokes_ind=0 #we want stokes I

rapplot.plot_data_stokes(image,min,max,stokes_ind,data_id,fig,axs[0][0],halfrange,mas,label="Stokes I",cmap="afmhot")

rapplot.plot_data_stokes(image,min,max,1,data_id,fig,axs[1][0],halfrange,mas,label="Stokes Q",cmap="RdBu")

rapplot.plot_data_stokes(image,min,max,2,data_id,fig,axs[0][1],halfrange,mas,label="Stokes U",cmap="RdBu")

rapplot.plot_data_stokes(image,min,max,3,data_id,fig,axs[1][1],halfrange,mas,label="Stokes V",cmap="RdBu")

fig.suptitle('t=%.01lf [days]'%(ind*10.*Tunit),fontsize=20)

axs[0][0].set_xlabel(r"x [mas]")
axs[0][0].set_ylabel(r"y [mas]")

axs[0][0].set_xlim(-0.1,0.1)
axs[0][0].set_ylim(-0.1,0.1)

axs[1][0].set_xlabel(r"x [mas]")
axs[1][0].set_ylabel(r"y [mas]")

axs[1][0].set_xlim(-0.1,0.1)
axs[1][0].set_ylim(-0.1,0.1)

axs[0][1].set_xlabel(r"x [mas]")
axs[0][1].set_ylabel(r"y [mas]")

axs[0][1].set_xlim(-0.1,0.1)
axs[0][1].set_ylim(-0.1,0.1)

axs[1][1].set_xlabel(r"x [mas]")
axs[1][1].set_ylabel(r"y [mas]")

axs[1][1].set_xlim(-0.1,0.1)
axs[1][1].set_ylim(-0.1,0.1)

plt.tight_layout()
Path("figures/").mkdir(parents=True, exist_ok=True)
print("figures/"+"img_%d.png"%ind)
plt.savefig("figures/img_%d.png"%ind, transparent=False,dpi=250)
plt.clf()

image.close()
