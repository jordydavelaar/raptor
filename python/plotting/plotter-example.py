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
M=6.5e9 * rapplot.MSUN
d=16800 * rapplot.KPC

rg = (rapplot.G*M/rapplot.SPEED_OF_LIGHT**2.)

mas = (rg/d)* rapplot.MAS_IN_DEG

Tunit =rg/rapplot.SPEED_OF_LIGHT

halfrange=200 #in rg

#read data
min,max,image = rapplot.read_data("output",ind,data_id)

#plot data
plt.figure(figsize=(6,5),dpi=250,facecolor='w')
fig, axs = plt.subplots(1,1,figsize=(6,5))

stokes_ind=0 #we want stokes I

rapplot.plot_data_stokes(image,min,max,stokes_ind,data_id,fig,axs,halfrange,mas,label="Stokes I",cmap="afmhot")

fig.suptitle('t=%.01lf [days]'%(ind*10.*Tunit),fontsize=20)

axs.set_xlabel(r"x [mas]")
axs.set_ylabel(r"y [mas]")

plt.tight_layout()
Path("figures/").mkdir(parents=True, exist_ok=True)
print("figures/"+"img_%d.png"%ind)
plt.savefig("figures/img_%d.png"%ind, transparent=False)
plt.clf()

image.close()
