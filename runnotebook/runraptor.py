import math
import numpy as np
import matplotlib.pyplot as plt
import emoji
from ipywidgets import Layout, interact, IntSlider,Label, HBox
import ipywidgets as widgets
import os
import sys
import h5py
from pathlib import Path
import rapplot

class pythonraptor():
	def __init__(self,image_width=None,image_height=None):
		self.inclination=160
		self.mass = 6.2e9
		self.distance = 3945
		self.Frequency = 23
		self.cam_x=40
		self.cam_y=40
		self.M_unit=2.739e25
		self.R_low=1
		self.R_high=1
		
		self.image_width = image_width
		self.image_height = image_height
		self.freqs_per_dec=1
		self.stepsize = 0.05
		self.max_level = 1
		
		if image_width==None:
			self.image_width = 200
		if image_height==None:
			self.image_height = 200



	def makeslider(self,imin=0,imax=180,minmass=1.e3,maxmass=1.e10,mindis=0.2,maxdis=4000,feqmin=10,feqmax=30,camxmin=15,camxmax=40,camymin=15,camymax=40):
		self.inclination = IntSlider(min=imin,max=imax,step=1) 
		inclination_box = HBox([Label('Inclination (Deg):'), self.inclination])
	
		self.mass = IntSlider(min=minmass,max=maxmass,step=1)
		mass_box = HBox([Label('Mass of Black hole in units of solar mass'), self.mass])
	
		self.distance = IntSlider(min=mindis,max=maxdis,step=1)
		distance_box = HBox([Label('Distance in $10^3$ parsec'), self.distance])
	
		self.Frequency = IntSlider(min=feqmin,max=feqmax,step=1)
		Frequency_box = HBox([Label('Frequency (x $10^{10}$) Hz:'), self.Frequency])
	
		self.cam_x = IntSlider(min=camxmin,max=camxmax,step=1)
		cam_x_box = HBox([Label('Camera Size (X):'), self.cam_x])

		self.cam_y = IntSlider(min=camymin,max=camymax,step=1)
		cam_y_box = HBox([Label('Camera Size (Y):'), self.cam_y])

		display(inclination_box)
		display(mass_box)
		display(distance_box)
		display(cam_x_box)
		display(cam_y_box)
		display(Frequency_box)
		return self.inclination,self.mass,self.distance,self.Frequency,self.cam_x,self.cam_y


	def writefile(self):

		freq_min = self.Frequency.value*1.e10
		#Writing the model.in file 
		f = open('model.in','r')
		text = f.readlines()
		f.close()
		text[0] = 'MBH\t\t(Msun)\t%.3e\n'%(self.mass.value)
		text[1] = 'DISTANCE\t\t(kpc)\t%.3e\n'%(self.distance.value)
		text[2] = 'M_UNIT\t\t(g)\t%.3e\n'%(self.M_unit)
		text[3] = 'R_LOW\t\t(-)\t%d\n'%(self.R_low)
		text[4] = 'R_HIGH\t\t(-)\t%d\n'%(self.R_high)
		text[5] = '\t\n'
		text[6] = 'INCLINATION\t(deg)\t\t%d\n'%(self.inclination.value)
		text[7] = 'IMG_WIDTH\t(pixels)\t%d\n'%(self.image_width)
		text[8] = 'IMG_HEIGHT\t(pixels)\t%d\n'%(self.image_height)
		text[9] = 'CAM_SIZE_X\t(Rg)\t\t%d\n'%(self.cam_x.value)
		text[10] = 'CAM_SIZE_Y\t(Rg)\t\t%d\n'%(self.cam_y.value)
		text[11] = 'FREQS_PER_DEC\t(-)\t\t%d\n'%(self.freqs_per_dec)
		text[12] = 'FREQ_MIN\t(Hz)\t\t%.3e\n'%(freq_min)
		text[13] = 'STEPSIZE\t(-)\t\t%.2e\n'%(self.stepsize)
		text[14] = 'MAX_LEVEL\t(-)\t\t%d\n'%(self.max_level)
		f = open('model.in','w')
		f.writelines(text)
		f.close()


	def plot(self,fileindex):
	
		M=self.mass.value * rapplot.MSUN
		d=self.distance.value * rapplot.KPC
		rg = (rapplot.G*M/rapplot.SPEED_OF_LIGHT**2.)
		mas = (rg/d)* rapplot.MAS_IN_DEG
		Tunit =rg/rapplot.SPEED_OF_LIGHT

		halfrange=int(self.cam_x.value/2)
		ind = rapplot.read_data_id('output',fileindex)
		min,max,image = rapplot.read_data('output',fileindex,ind)
		plt.figure(figsize=(6,5),dpi=500,facecolor='w')
		fig, ax = plt.subplots(1,1,figsize=(6,5))
		p = rapplot.plot_data_stokes(image,min,max,0,ind,fig,ax,halfrange,mas,label="Stokes I",cmap="afmhot")

		ax.set_xlabel(r"x [mas]")
		ax.set_ylabel(r"y [mas]")
		plt.title('t=%.01lf [days]'%(ind*10.*Tunit),fontsize=20)

		plt.tight_layout()
		plt.show()
		return 0
