'''
Set parameters in 'MultiNestRAPTORinput.txt'
Run this file with command: python3 Multinest.py Run-name

To run this code, one has to have downloaded PyMultiNest, eht-imaging

Needed files:
	Multinest.py
	MultiNestRAPTORinput.txt
	Convergenceplotter.py

Maybe needed files:
	Observational_Spectrum.txt
	Observational_Coreshift.txt
        Observational_Image.txt

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

Developped by:
Joost de Kleuver, Renze Oosterhuis

De oplossing van de formule? Moleculen zonder scrupules
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
import ehtim as eh

# Define global variables
MBH_var, M_UNIT_var, Rhigh_var, i_var, Spectrum, Coreshift, Image = [0,0,0,0,0,0,0]
num_Spectrum, num_IR, num_Coreshift, num_indices = [0,0,0,0]
MBH, M_UNIT, Rhigh, i, data_number = [0,0,0,0,0]
loglike_evaluation = 0

# Define functions
def Initialize_imgrendererc():
	global num_Spectrum, num_Coreshift
        if Spectrum: # Radio frequencies
                FREQ_Spectrum = np.loadtxt('Observational_Spectrum.txt').transpose()[0]
        	num_Spectrum = len(FREQ_Spectrum)
        if Coreshift: # Coreshift frequencies
                FREQ_Coreshift = np.loadtxt('Observational_Coreshift.txt').transpose()[0]
#                FREQ_Coreshift = np.divide(29979245800,np.loadtxt('Observations_Coreshift.txt').transpose()[0])
		num_Coreshift = len(FREQ_Coreshift)
        FREQ = np.concatenate([FREQ_Spectrum,FREQ_Coreshift])
	global FREQ
	f = open('frequencies.txt', 'w')
	for i in FREQ:
		f.write(str(i + "\n")
	f.close()
			
def Initialize_sizeflux_v4c(IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, source_dist):
        if Coreshift: # Coreshift frequencies
                FREQ_Coreshift = np.loadtxt('Observations_Coreshift.txt').transpose()[0]
#                FREQ_Coreshift = np.divide(29979245800,np.loadtxt('Observations_Coreshift.txt').transpose()[0])
		f = open('sizeflux_v4_5.c','r')
        	text = f.readlines()
        	f.close()
		text[11] = '\tint dump1,dump2,img_size1=%d,img_size2=%d;\n'%(IMG_WIDTH, IMG_HEIGHT)
		text[23] = '\tdouble freqlist[%d] = {%s};\n'%(len(FREQ_Coreshift),', '.join(map(str, FREQ_Coreshift)))
		text[40] = '\tdouble d=%.15e;//5.061e25;\n'%(source_dist)
        	text[43] = '\tdouble CAM_SIZE_X = %f;\n'%(CAM_SIZE_X)
		f = open('sizeflux_v4_5.c','w')
        	f.writelines(text)
        	f.close()
		os.system('gcc -o sizeflux_v4_5 sizeflux_v4_5.c -lm')
			
def Initialize_modelin(params_modelin):
	# Set parameters in 'model.in'
	MBH, M_UNIT, Rlow, Rhigh, i, IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, CAM_SIZE_Y, FREQS_PER_DEC, FREQ_MIN, STEPSIZE, MAX_LEVEL = params_modelin
	f = open('model.in','r')
	text = f.readlines()
	f.close()
	text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
	text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
	text[2] = 'Rlow\t\t(-)\t\t%d\n'%(Rlow)
	text[3] = 'Rhigh\t\t(-)\t\t%d\n'%(Rhigh)

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

def Initialize_parametersh(source_dist):
    	# Set parameters in 'parameters.h'
        f = open('parameters.h','r')
        text = f.readlines()
        f.close()
	num_frequencies = num_Spectrum + num_Coreshift
        text[37] = '#define num_indices %d\n'%(num_frequencies)
	text[42] = '#define FREQS (FREQFILE) \n'
	text[142] = '#define source_dist    (%.15e) // Distance to Centaurus A (cm)\n'%(source_dist)
        f = open('parameters.h','w')
        f.writelines(text)
        f.close()
	os.system('make clean')
	os.system('make all')
			
def Set_parameters():
    	# Load parameters from 'MultiNestRAPTORinput.txt'
	_, params = np.loadtxt('MultiNestRAPTORinput.txt', dtype=str, delimiter='=', unpack=True)
	params = params[0:22].astype(float)
       	# Set global variables
	global MBH_var, M_UNIT_var, Rhigh_var, i_var, Spectrum, Image, Coreshift
	MBH_var, M_UNIT_var, Rhigh_var, i_var, Spectrum, Image, Coreshift = params[0:7].astype(int)
	global MBH, M_UNIT, Rhigh, i, data_number, source_dist
	if not MBH_var:
		MBH = params[7]
	if not M_UNIT_var:
		M_UNIT = params[8]
	if not Rhigh_var:
		Rhigh = params[9]
	if not i_var:
		i = params[11]
	source_dist = params[20]
	data_number = params[22]
    	# Initialize input files RAPTOR
	Initialize_modelin(params[7:20])

def Set_modelin(MBH, M_UNIT,Rhigh,i):
        # Set MBH, M_UNIT, Rhigh and i in 'model.in'
        f = open('model.in','r')
        text = f.readlines()
        f.close()
        text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
        text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
        text[3] = 'Rhigh\t\t(-)\t\t%d\n'%(Rhigh)
        text[5] = 'INCLINATION\t(deg)\t%d\n'%(i)
        f = open('model.in','w')
        f.writelines(text)
        f.close()

def RAPTOR(MBH, M_UNIT, Rhigh, i, data_number):
        # Set parameter values for RAPTOR:
        Set_modelin(MBH, M_UNIT, Rhigh,i)
    	# Run RAPTOR
        os.system('./RAPTOR model.in ../data1000.dat %d'%(data_number))
    	# Load flux values
	Chisquare_Im = 0
        if Image:
		obs = eh.obsdata.load_uvfits(‘obsname.uvfits’)
		im = eh.image.load_txt(‘output/uniform_img_%d_%d’%(FREQ[0],data_number))
		im.rf = obs.rf
		im.ra = obs.ra
		im.dec = obs.dec
		#get reduced chi squares
		chi_cphase  = obs.chisq(im, dtype= ‘cphase’’)
		chi_logcamp = obs.chisq(im, dtype= ‘logcamp’)
		Chisquare_Im = chi_cphase + chi_logcamp
        Freq, Flux = np.loadtxt('output/spectrum_%d_%d.00.dat'%(data_number,np.floor(i)), dtype=str, unpack=True)
        Freq = Freq.astype(float)
        Flux = Flux.astype(float)
    	# Write output to 'Run-name_flux.txt'
        f = open("%s_Flux.txt"%(sys.argv[1]),'a+')
        f.write(str(Freq)+' '+str(Flux)+'\n')
        # f.write(' '.join(map(str,Flux))+'\n')
        f.close()
    	# Return flux values
        # 
        return Flux, Chisquare_Im

def myprior(cube, ndim=MBH_var+M_UNIT_var+Rhigh_var+i_var, nparams=MBH_var+M_UNIT_var+Rhigh_var+i_var):
    	if MBH_var: # MBH has a log uniform distribution between 1.0e40 and 1.0e42 g
        	cube[MBH_var-1] = 10**(2.0*cube[MBH_var-1] + 39)
    	if M_UNIT_var: # M_unit has a log uniform distribution between 1e22 and 1e24 g
        	cube[MBH_var+M_UNIT_var-1] = 10**(2*cube[MBH_var+M_UNIT_var-1] + 22)
    	if Rhigh_var: # Rhigh has a uniform distribution between 1 and 100
        	cube[MBH_var+M_UNIT_var+Rhigh_var-1] = 99*cube[MBH_var+M_UNIT_var+Rhigh_var-1] + 1
    	if i_var: # i has a uniform distribution between 12 and 45
        	cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] = 33*cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] + 12
    	return cube

def myloglike(cube, ndim=MBH_var+M_UNIT_var+Rhigh_var+i_var, nparams=MBH_var+M_UNIT_var+Rhigh_var+i_var):
        global loglike_evaluation
        loglike_evaluation += 1
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
        #Calculate Chi^2
        Flux, Chisquare_Image = RAPTOR(MBH, M_UNIT, Rhigh, i, data_number)
        if Spectrum:
                _, Measurements_Spectrum, Sigma_Spectrum = np.loadtxt('Observational_Spectrum.txt', delimiter = " ", unpack = True)
        	global num_Spectrum
        	for freqnum in range(num_Spectrum):
#            		Chisquare_Radio += (Measurements_Radio[freqnum] - np.log10(Flux[freqnum]))**2  /(2*Sigma_Radio[freqnum]**2)
        		Chisquare_Spectrum += (Measurements_Spetrum[freqnum] - Flux[freqnum])**2  /(num_Spectrum*2*Sigma_Spectrum[freqnum]**2)
#			Chisquare_Radio += (np.log10(Measurements_Radio[freqnum]) - np.log10(Flux[freqnum]))**2 /(2*(np.log10(Measurements_Radio[freqnum])-np.log10(Measurements_Radio[freqnum] - Sigma_Radio[freqnum]))**2) #       	
#			Chisquare_Radio +=(Measurements_Radio[freqnum] - Flux[freqnum])**2
        # Calculate total chisquared value
    	if Coreshift:
#		_, Measurements_Major, Sigma_Major, Measurements_Minor, Sigma_Minor = np.loadtxt('Observations_Coreshift.txt').transpose()
		_, Measurements_Major, Err_up_Major, Err_down_Major = np.loadtxt('Observations_Coreshift.txt').transpose()
        	os.system('./sizeflux_v4_5 1 %d %.2f %f %d'%(num_Coreshift, i, MBH, data_number))
            	Wavelength, Major, Minor, _, _, _, _ = np.loadtxt('output/lambda_th_%d_%.2f.dat'%(data_number,i), dtype=float, delimiter=' ', unpack=True)
		for Wavelengthnum in range(len(Wavelength)):
                	Chisquare_Coreshift += (Measurements_Major[Wavelengthnum] - Major[Wavelengthnum])**2  /(num_Coreshift*2*Err_up_Major[Wavelengthnum]**2)
#			Chisquare_Coreshift += (Measurements_Minor[Wavelengthnum] - Minor[Wavelengthnum])**2  /(2*Sigma_Minor[Wavelengthnum]**2)
#			Chisquare_Coreshift += (np.log10(Measurements_Major[Wavelengthnum]) - np.log10(Major[Wavelengthnum]))**2 /(2*(np.log10(Measurements_Coreshift[Wavelenghtnum])-np.log10(Measurements_Coreshift[Wavelengthnum] - Err_down_Major[Wavelengthnum]))**2)
		# Write output to 'Run-name_Major.txt'
		
	Chisquare_Total = Chisquare_Spectrum + Chisquare_Image + Chisquare_Coreshift
    	# Write Chisquared values to 'Run-name_Chisquared.txt'
        f = open("%s_Chisquared.txt"%(sys.argv[1]),'a+')
        f.write("%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n"%(MBH, M_UNIT, Rhigh, i, Chisquare_Total))
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
        os.system('python3 Convergenceplotter.py %s'%(sys.argv[1]))
    	# Return loglikelihood
        return - Chisquare_Total

## Main code
Set_parameters()

parameters = np.array(["MBH", "M_UNIT", "Rhigh", "i"])[np.array([MBH_var,M_UNIT_var,Rhigh_var,i_var]).astype(bool)]
print(parameters)
n_params = len(parameters)
prefix = sys.argv[1]

result = pymultinest.run(LogLikelihood=myloglike, Prior=myprior, max_iter = 0,n_dims=n_params, importance_nested_sampling=False, multimodal=True, evidence_tolerance=0.5, sampling_efficiency=0.8, null_log_evidence=-1e90,outputfiles_basename = prefix, max_modes=5, mode_tolerance=-1e90, resume=False, verbose=True, write_output=True, n_live_points=400)

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

