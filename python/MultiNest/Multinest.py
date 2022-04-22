'''
Set parameters in 'MultiNestRAPTORinput.txt'
Check if there are enough lines for the frequencies in file 'img_renderer.c'
Run this file with command: python Multinest.py Run-name

To run this code, one has to have downloaded PyMultiNest

Needed files:
        Multinest.py
	MultiNestRAPTORinput.txt
	Convergenceplotter.py
	
	spectrum.txt

Not used in this code, but needed to make plots of the best fits afterwards
	plotterSgrAalles.py

Output consists of:
	Run-name_Flux.txt
	Run-name_Major.txt
	Run-name_Minor.txt
	Run-name_Chisquared.txt
	Run-name_Convergence.png
	Run-name_livepoints.txt
	Run-name_modes.txt
'''
#====================================================================================================
# Import packages
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pymultinest
from pymultinest.solve import solve
from pymultinest.analyse import Analyzer

# Define global variables
MBH_var, M_UNIT_var, Rhigh_var, i_var, Radio, IR, Coreshift = [0,0,0,0,0,0,0]
num_Radio, num_IR, num_Coreshift, num_indices = [0,0,0,0]
MBH, M_UNIT, Rhigh, i, data_number = [0,0,0,0,0]
loglike_evaluation = 0
data_number = 0
# Define functions
def Initialize_modelin(params_modelin):
	# Set parameters in 'model.in'
	MBH, M_UNIT, Rhigh, Rlow, i, IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, CAM_SIZE_Y, FREQS_PER_DEC, FREQ_MIN, STEPSIZE, MAX_LEVEL = params_modelin
	f = open('model.in','r')
	text = f.readlines()
	f.close()
	text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
	text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
	text[2] = 'Rhigh\t\t(-)\t\t%d\n'%(Rhigh)
	text[3] = 'Rlow\t\t(-)\t\t%d\n'%(Rlow)
	text[5] = 'INCLINATION\t(deg)\t%d\n'%(i)
	text[6] = 'IMG_WIDTH\t(pixels)\t%d\n'%(IMG_WIDTH)
	text[7] = 'IMG_HEIGHT\t(pixels)\t%d\n'%(IMG_HEIGHT)
	text[8] = 'CAM_SIZE_X\t(Rg)\t\t%.15e\n'%(CAM_SIZE_X)
	text[9] = 'CAM_SIZE_Y\t(Rg)\t\t%.15e\n'%(CAM_SIZE_Y)
	text[10] = 'FREQS_PER_DEC\t(-)\t\t%d\n'%(FREQS_PER_DEC)
	text[11] = 'FREQ_MIN\t(Hz)\t\t%.15e\n'%(FREQ_MIN)
	text[12] = 'STEPSIZE\t(-)\t\t%.15e\n'%(STEPSIZE)
	text[13] = 'MAX_LEVEL\t(-)\t\t%d\n'%(MAX_LEVEL)
	f = open('model.in','w')
	f.writelines(text)
	f.close()

def Set_parameters():
    	# Load parameters from 'MultiNestRAPTORinput.txt'
	_, params = np.loadtxt('MultiNestRAPTORinput.txt', dtype=str, delimiter='=', unpack=True)
	params = params.astype(float)
       	# Set global variables
	global MBH_var, M_UNIT_var, Rhigh_var, i_var, Radio, IR, Coreshift
	MBH_var, M_UNIT_var, Rhigh_var, i_var = params[0:4].astype(int)
	global MBH, M_UNIT, Rhigh, i, data_number
	if not MBH_var:
		MBH = params[4]
	if not M_UNIT_var:
		M_UNIT = params[5]
	if not Rhigh_var:
		Rhigh = params[6]
	if not i_var:
		i = params[8]
	#data_number = params[26]
	data_number = 0
    	# Initialize input files RAPTOR
	Initialize_modelin(params[4:16])

def Set_modelin(MBH, M_UNIT,Rhigh,i):
        # Set MBH, M_UNIT, Rhigh and i in 'model.in'
        f = open('model.in','r')
        text = f.readlines()
        f.close()
        text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
        text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
        text[2] = 'Rhigh\t\t(-)\t\t%d\n'%(Rhigh)
        text[5] = 'INCLINATION\t(deg)\t%d\n'%(i)
        f = open('model.in','w')
        f.writelines(text)
        f.close()

def RAPTOR(MBH, M_UNIT, Rhigh, i, data_number):
        # Set parameter values for RAPTOR:
        Set_modelin(MBH, M_UNIT, Rhigh,i)
    	# Run RAPTOR
        os.system('./RAPTOR model.in ../harm3d.txt %d'%(data_number))
    	# Load flux values
        Freq, Flux = np.loadtxt('output/spectrum_%d_%d.00.dat'%(data_number,np.floor(i)), dtype=str, unpack=True)
        Freq = Freq.astype(float)
        Flux = Flux.astype(float)
    	# Write output to 'Run-name_flux.txt'
        f = open("%s_Flux.txt"%(sys.argv[1]),'a+')
        f.write(str(Freq)+' '+str(Flux)+'\n')
        f.close()
    	# Return flux values
        return Flux

def myprior(cube, ndim=MBH_var+M_UNIT_var+Rhigh_var+i_var, nparams=MBH_var+M_UNIT_var+Rhigh_var+i_var):
    	if MBH_var: # MBH has a log uniform distribution between 1.0e39 and 1.0e41 g
        	cube[MBH_var-1] = 10**(2.0*cube[MBH_var-1] + 39)
    	if M_UNIT_var: # M_unit has a log uniform distribution between 1e19 and 1e24 g
        	cube[MBH_var+M_UNIT_var-1] = 10**(5*cube[MBH_var+M_UNIT_var-1] + 19)
    	if Rhigh_var: # Rhigh has a log uniform distribution between 10 and 100
        	cube[MBH_var+M_UNIT_var+Rhigh_var-1] = 10**(cube[MBH_var+M_UNIT_var+Rhigh_var-1] + 1)
    	if i_var: # i has a uniform distribution between 12 and 45
        	cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] = 33*cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] + 12
    	return cube

def myloglike(cube, ndim=MBH_var+M_UNIT_var+Rhigh_var+i_var, nparams=MBH_var+M_UNIT_var+Rhigh_var+i_var):
        global loglike_evaluation
        global data_number
        loglike_evaluation += 1
        data_number += 1
	# Get flux values
        global MBH, M_UNIT, Rhigh, i
        if MBH_var:
                MBH = cube[MBH_var-1]
        if M_UNIT_var:
                M_UNIT = cube[MBH_var+M_UNIT_var-1]
        if Rhigh_var:
                Rhigh = cube[MBH_var+M_UNIT_var+Rhigh_var-1]
        if i_var:
                i = cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1]
        Flux = RAPTOR(MBH, M_UNIT, Rhigh, i, data_number)
    	# Calculate chisquared values for radio, IR and Coreshift
        _, Measurements_Radio = np.loadtxt('spectrum.txt').transpose()
        global num_Radio
        num_Radio = 8
        Chisquare_Radio = 0
        for freqnum in range(num_Radio):
#            	Chisquare_Radio += (Measurements_Radio[freqnum] - np.log10(Flux[freqnum]))**2  /(2*Sigma_Radio[freqnum]**2)
                Chisquare_Radio += (Measurements_Radio[freqnum] - Flux[freqnum])**2
#		Chisquare_Radio += (np.log10(Measurements_Radio[freqnum]) - np.log10(Flux[freqnum]))**2 /(2*(np.log10(Measurements_Radio[freqnum])-np.log10(Measurements_Radio[freqnum] - Sigma_Radio[freqnum]))**2)
    	# Write Chisquared values to 'Run-name_Chisquared.txt'
        f = open("%s_Chisquared.txt"%(sys.argv[1]),'a+')
        f.write("%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n"%(MBH, M_UNIT, Rhigh, i, Chisquare_Radio))
        f.close()
	# Write all live point info to 'Run-name_livepoints.txt'
        livepoints = np.loadtxt('%sphys_live.points'%(sys.argv[1]), dtype=str)
        f = open("%s_livepoints.txt"%(sys.argv[1]),'a+')
        for i in range(len(livepoints)):
                livepoint_i = livepoints[i]
                livepoint_i = np.append(livepoint_i, np.array([loglike_evaluation]))
                f.write(str(livepoint_i) + '\n')
        f.close()
	# Update Convergence plot
        #os.system('python Convergenceplotter.py %s'%(sys.argv[1]))
    	# Return loglikelihood
        return Chisquare_Radio

## Main code
Set_parameters()

parameters = np.array(["MBH", "M_UNIT", "Rhigh", "i"])[np.array([MBH_var,M_UNIT_var,Rhigh_var,i_var]).astype(bool)]
print(parameters)
n_params = len(parameters)
prefix = sys.argv[1]

result = pymultinest.run(LogLikelihood=myloglike, Prior=myprior,n_dims=n_params, importance_nested_sampling=False, multimodal=True, evidence_tolerance=0.5, sampling_efficiency=0.8, null_log_evidence=-1e90,outputfiles_basename = prefix, max_iter =5,  max_modes=5, mode_tolerance=-1e90, resume=False, verbose=True, write_output=True, n_live_points=100)
print("gelukt")

# Write for each mode the mean values of the parameters with the standard deviation and the log_evidence with its error to prefix_modes.txt
a = Analyzer(n_params, outputfiles_basename="%s"%(prefix))
stats = a.get_stats()
f = open("%s_modes.txt"%(prefix),'a+')
for mode in stats[u'modes']:
	mean = mode[u'mean']
	sigma = mode[u'sigma']
	log_evidence = mode[u'local log-evidence']
	log_evidence_error = mode[u'local log-evidence error']
	string = ''
	for param in range(n_params):
		string += '%.15e\t%.15e\t'%(mean[param],sigma[param])
	string += '%.15e\t%.15e\n'%(log_evidence,log_evidence_error)
	f.write(string)
f.close()


# Best fit can be found by calculating the chisquared for each resulting parameter estimate and choosing the lowest one
# SED fit and core size plot can be made using plotterSgrAalles.py

