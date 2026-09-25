import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt
import h5py

G = 6.674e-8
MSUN = 1.989e33
SPEED_OF_LIGHT = 2.998e10

KPC = 3.086e21
SEC_IN_DAY = 86400.
MAS_IN_DEG = 206264.806*1000.

def read_data_id(folder,ind):
    file_name = folder+'/img_data_%d.h5'%ind
    print("Reading keys from: ", file_name)
    images = h5py.File(file_name,'r')
    keys = [key for key in images.keys()]
    print(keys)
    images.close()
    return keys

def read_data(folder,ind,data_id):

    file_name = folder+'/img_data_%d.h5'%ind

    print("Reading in: ", file_name)

    images = h5py.File(file_name,'r')

    min = [-100.,-100.,-100.,-100.]
    max = [100.,100.,100,100.]
    print(len(data_id))
    print(data_id)
    for j in range(0,4):
        for i in range(0,len(images[data_id[j]])):
            current=np.max(images[data_id[j]][i])
            max[j]=np.maximum(max[j],np.max(images[data_id[j]][i]))
            min[j]=np.minimum(min[j],np.min(images[data_id[j]][i]))


    return min,max,images

def plot_data_tau(image,data_id,ind,fig,ax,halfrange=40,mas=1,label="Stokes",cmap="CMRmap",vmin=-8,vmax=2):

    for i in range(0,len(image[data_id[ind]])):
        pixels=int(np.sqrt(len(image[data_id[ind]][i])))
        array=((np.reshape(image[data_id[ind]][i],(pixels,pixels))))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas

        figure=ax.pcolormesh(alpha,beta,np.log10(array+1e-10),vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')
        ax.set_aspect('equal')
        
    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_stokes(image,min,max,stokes_ind,data_id,fig,ax,halfrange=40,mas=1,label="Stokes",cmap="afmhot"):

    for i in range(0,len(image[data_id[stokes_ind]])):
        pixels=int(np.sqrt(len(image[data_id[stokes_ind]][i])))
        array=((np.reshape(image[data_id[stokes_ind]][i],(pixels,pixels))))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas
        ax.set_aspect('equal')

        if(stokes_ind==0):
            figure=ax.pcolormesh(alpha,beta,(array/max[stokes_ind])**0.5,vmin=0,vmax=1,cmap=cmap,shading='auto')
        else:
            figure=ax.pcolormesh(alpha,beta,(array/max[stokes_ind]),vmin=-1,vmax=1,cmap=cmap,shading='auto')

    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_norder(image,fig,ax,halfrange=40,mas=1,label="n (winding number)",cmap="RdBu_r",vmin=-3,vmax=3):

    for i in range(0,len(image['norder'])):
        pixels=int(np.sqrt(len(image['norder'][i])))
        array=((np.reshape(image['norder'][i],(pixels,pixels))))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas
        ax.set_aspect('equal')

        figure=ax.pcolormesh(alpha,beta,array,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')

    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_norder_mino(image,fig,ax,halfrange=40,mas=1,label="n (Mino-time half-orbit count)",cmap="RdBu_r",vmin=-3,vmax=3):
    # Gralla & Lupsasca 2020 (arXiv:1910.12873) / Himwich, Johnson, Lupsasca &
    # Strominger 2020 (arXiv:2001.08750) convention. Sentinel value -999
    # marks vortical geodesics (eta<=0) where this is undefined; those
    # pixels are masked out rather than plotted.
    for i in range(0,len(image['norder_mino'])):
        pixels=int(np.sqrt(len(image['norder_mino'][i])))
        array=((np.reshape(image['norder_mino'][i],(pixels,pixels))))
        array=np.ma.masked_less(array,-100)
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas
        ax.set_aspect('equal')

        figure=ax.pcolormesh(alpha,beta,array,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')

    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_ncross(image,fig,ax,halfrange=40,mas=1,label="equatorial crossings",cmap="tab10",vmin=0,vmax=9):
    for i in range(0,len(image['ncross'])):
        pixels=int(np.sqrt(len(image['ncross'][i])))
        array=((np.reshape(image['ncross'][i],(pixels,pixels))))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas
        ax.set_aspect('equal')

        figure=ax.pcolormesh(alpha,beta,array,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')

    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_dchi_grav(image,fig,ax,halfrange=40,mas=1,label=r"$\Delta\chi_{\rm grav}$ [deg]",cmap="twilight_shifted",vmin=-90,vmax=90):
    # EVPA rotation from the sky-projected spin axis between a source at
    # infinity and the camera (dchi_grav, radians in the file). NaN for rays
    # that are captured or stopped early; pcolormesh leaves those blank.
    for i in range(0,len(image['dchi_grav'])):
        pixels=int(np.sqrt(len(image['dchi_grav'][i])))
        array=np.degrees(np.reshape(image['dchi_grav'][i],(pixels,pixels)))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas
        ax.set_aspect('equal')

        figure=ax.pcolormesh(alpha,beta,array,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')

    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def overlay_norder_contours(image,ax,levels=(0.75,1.25,2.0,3.0),mas=1):

    for i in range(0,len(image['norder'])):
        pixels=int(np.sqrt(len(image['norder'][i])))
        array=((np.reshape(image['norder'][i],(pixels,pixels))))
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas

        ax.contour(alpha,beta,np.abs(array),levels=levels,colors='cyan',linewidths=0.5)

def overlay_norder_mino_contours(image,ax,levels=(1,2,3,4,5),mas=1,colors='cyan',label=True):
    # norder_mino is already a half-orbit count (n=1,2,3,... at the
    # turning-point boundaries between subrings), so contour directly at
    # integer levels rather than at offset thresholds like the winding
    # number. Sentinel -999 (vortical geodesics) is masked out. The shadow
    # is generally offset/asymmetric for spinning Kerr holes at nonzero
    # inclination, so contours are not concentric circles -- label=True
    # annotates each line with its n value directly rather than relying on
    # radius/ordering by eye.
    for i in range(0,len(image['norder_mino'])):
        pixels=int(np.sqrt(len(image['norder_mino'][i])))
        array=((np.reshape(image['norder_mino'][i],(pixels,pixels))))
        array=np.ma.masked_less(array,-100)
        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((-np.reshape(image['beta'][i],(pixels,pixels))))*mas

        cs=ax.contour(alpha,beta,array,levels=levels,colors=colors,linewidths=0.5)
        if label:
            ax.clabel(cs,levels,inline=True,fontsize=7,fmt='n=%d')

def plot_data_polfrac(image,max,data_id,fig,ax,halfrange=10,mas=1,label="|m|",cmap="afmhot"):

    for i in range(0,len(image[data_id[0]])):
        pixels=int(np.sqrt(len(image[data_id[0]][i])))

        array_I=((np.reshape(image[data_id[0]][i],(pixels,pixels))))
        array_Q=((np.reshape(image[data_id[1]][i],(pixels,pixels))))
        array_U=((np.reshape(image[data_id[2]][i],(pixels,pixels))))

        array = np.sqrt(array_Q**2.+array_U**2)/array_I
        array[array_I/max[0]<1e-7]=0

        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((np.reshape(-image['beta'][i],(pixels,pixels))))*mas

        figure=ax.pcolormesh(alpha,beta,(array),vmin=0,vmax=1,cmap=cmap,shading='auto')


    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)

def plot_data_RM(image,min,max,ind_1,ind_2,data_id,lam1,lam2,fig,ax,halfrange=10,mas=1,label="RM",cmap="RdBu"):

    lam1*=1e-3
    lam2*=1e-3
    for i in range(0,len(image[data_id[0]])):
        pixels=int(np.sqrt(len(image[data_id[0]][i])))

       	array_I_1=((np.reshape(image[data_id[ind_1]][i],(pixels,pixels))))
        array_Q_1=((np.reshape(image[data_id[ind_1 + 20]][i],(pixels,pixels))))
        array_U_1=((np.reshape(image[data_id[ind_1+ 40]][i],(pixels,pixels))))

        array_I_2=((np.reshape(image[data_id[ind_2]][i],(pixels,pixels))))
        array_Q_2=((np.reshape(image[data_id[ind_2+20]][i],(pixels,pixels))))
        array_U_2=((np.reshape(image[data_id[ind_2+40]][i],(pixels,pixels))))

        LP_1 = np.sqrt(array_Q_1**2. + array_U_1**2)
        LP_2 = np.sqrt(array_Q_2**2. + array_U_2**2)

        EVPA_1 = 0.5*np.angle(array_Q_1+1j*array_U_1)

        EVPA_2 = 0.5*np.angle(array_Q_2+1j*array_U_2)
        EVPA_2[EVPA_1==0.0]=0
        EVPA_1[EVPA_2==0.0]=0

        array= EVPA_2
        for j in range(0,pixels):
           for k in range(0,pixels):
             array[j,k] = EVPA_2[j,k]-EVPA_1[j,k]
             if(np.abs(array[j,k])>0.75*np.pi):
                 array[j,k]=array[j,k]-np.sign(array[j,k])*np.pi

        RM = (array) #/(lam2**2-lam1**2)

        RM[array_I_1/max[0]<1e-6]=0
        RM[array_I_2/max[0]<1e-6]=0

        alpha=((np.reshape(image['alpha'][i],(pixels,pixels))))*mas
        beta=((np.reshape(-image['beta'][i],(pixels,pixels))))*mas

        figure=ax.pcolormesh(alpha,beta,np.sign(RM)*np.abs(RM)**0.25,vmin=-1,vmax=1,cmap=cmap,shading='auto')

        ax.set_xlim(-halfrange*mas,halfrange*mas)
        ax.set_ylim(-halfrange*mas,halfrange*mas)


    fig.colorbar(figure,label=label,ax=ax)

    ax.set_xlim(-halfrange*mas,halfrange*mas)
    ax.set_ylim(-halfrange*mas,halfrange*mas)
