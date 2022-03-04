#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt
import h5py

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}

matplotlib.rc('font', **font)


offset=0

def plot_data(folder,tstart,tend):
    plt.figure(figsize=(30,30),dpi=250,facecolor='w')
    fig, axs = plt.subplots(2,2,figsize=(30,30))

#    data_id = 'I3.425455e+15'
    data_I_id =  'I2.300000e+11'
    data_Q_id =  'Q2.300000e+11'
    data_U_id =  'U2.300000e+11'
    data_V_id =  'V2.300000e+11'

    for file in range(tstart,tend):
        t=(file+offset)%180
        print(file,t)
        file_name = folder+'/img_data_%d.h5'%t
        print(file_name)
        images = h5py.File(file_name,'r')
        print(images.keys())
        max_I=-100
        min_I=100
        max_Q=-100
        min_Q=100
        max_U=-100
        min_U=100
        max_V=-100
        min_V=100

        for i in range(0,len(images[data_I_id])):
            current_I=np.max(images[data_I_id][i])
            if(max_I<current_I):
                max_I=current_I
            current_I=np.min(images[data_I_id][i])
            if(min_I>current_I):
                min_I=current_I

            current_Q=np.max(images[data_Q_id][i])
            if(max_Q<current_Q):
                max_Q=current_Q
            current_Q=np.min(images[data_Q_id][i])
            if(min_Q>current_Q):
                min_Q=current_Q

            current_U=np.max(images[data_U_id][i])
            if(max_U<current_U):
                max_U=current_U
            current_U=np.min(images[data_U_id][i])
            if(min_U>current_U):
                min_U=current_U

            current_V=np.max(images[data_V_id][i])
            if(max_V<current_V):
                max_V=current_V
            current_V=np.min(images[data_V_id][i])
            if(min_V>current_V):
                min_V=current_V

        for i in range(0,len(images[data_I_id])):
            pixels=int(np.sqrt(len(images[data_I_id][i])))
            array_I=((np.reshape(images[data_I_id][i],(pixels,pixels))))
            array_Q=((np.reshape(images[data_Q_id][i],(pixels,pixels))))
            array_U=((np.reshape(images[data_U_id][i],(pixels,pixels))))
            array_V=((np.reshape(images[data_V_id][i],(pixels,pixels))))
            alpha=((np.reshape(images['alpha'][i],(pixels,pixels))))
            beta=((np.reshape(images['beta'][i],(pixels,pixels))))

            figure_I=axs[0][0].pcolormesh(alpha,beta,np.sqrt(array_I/max_I),vmin=0,vmax=1,cmap='afmhot')
            figure_Q=axs[1][0].pcolormesh(alpha,beta,(array_Q/max_Q),vmin=-1,vmax=1,cmap='RdBu',shading='auto')
            figure_U=axs[0][1].pcolormesh(alpha,beta,(array_U/max_U),vmin=-1,vmax=1,cmap='RdBu',shading='auto')
            figure_V=axs[1][1].pcolormesh(alpha,beta,(array_V/(max_V)),vmin=-0.00001,vmax=0.00001,cmap='RdBu',shading='auto')

        #plt.colorbar(figure_I,label="I")
        #plt.colorbar(figure_Q,label="I")
        #plt.colorbar(figure_U,label="I")
        #plt.colorbar(figure_V,label="I")

        plt.xlabel(r"x [r$_g$]")
        plt.ylabel(r"y [r$_g$]")

        #plt.axis('equal')
        #plt.tight_layout()
        Path("figures/"+folder).mkdir(parents=True, exist_ok=True)
        print("figures/"+folder+"/img_"+folder+"_%d.png"%file)
        plt.savefig("figures/"+folder+"/img_"+folder+"_%d.png"%file, transparent=False)
        plt.clf()
        images.close()

folder="model"
print(folder)
plot_data(folder,0,1)
