'''
Set parameters in 'MultiNestRAPTORinput.txt'
Check if there are enough lines for the frequencies in file 'img_renderer.c'
Run this file with command: python MultiNestRAPTOR_SgrA.py Run-name

To run this code, one has to have downloaded PyMultiNest

Needed files:
        MultiNestRAPTOR_SgrA.py
        sizeflux_v4_5.c   (and compile it with name sizeflux_v4_5)
	MultiNestRAPTORinput.txt
	Convergenceplotter.py

Maybe needed files:
	Observations_Radio.txt
	Observations_IR.txt
	Observations_Coreshift.txt

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

# Define functions
def Initialize_imgrendererc():
	global num_Radio, num_IR, num_Coreshift
        if Radio: # Radio frequencies
                FREQ_Radio = np.loadtxt('Observations_Radio.txt').transpose()[0]
        	num_Radio = len(FREQ_Radio)
	if IR: # IR frequencies
                FREQ_IR = np.loadtxt('Observations_IR.txt').transpose()[0]
		num_IR = len(FREQ_IR)
        if Coreshift: # Coreshift frequencies
#                FREQ_Coreshift = np.loadtxt('Observations_Coreshift.txt').transpose()[0]
                FREQ_Coreshift = np.divide(29979245800,np.loadtxt('Observations_Coreshift.txt').transpose()[0])
		num_Coreshift = len(FREQ_Coreshift)
        FREQ = np.concatenate([FREQ_Radio,FREQ_IR,FREQ_Coreshift])
        f = open('img_renderer.c','r')
        text = f.readlines()
        f.close()
        for freq_num in range(len(FREQ)):
                rule = 198 + freq_num
                text[rule-1] = '\tfrequencies[%d] = %.6e;\n'%(freq_num, FREQ[freq_num])
        for rule in range(198+len(FREQ), 250):
                text[rule-1] = '//\tfrequencies[99] = 0;\n'
        f = open('img_renderer.c','w')
        f.writelines(text)
        f.close()

def Initialize_sizeflux_v4c(IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, source_dist):
        if Coreshift: # Coreshift frequencies
#                FREQ_Coreshift = np.loadtxt('Observations_Coreshift.txt').transpose()[0]
                FREQ_Coreshift = np.divide(29979245800,np.loadtxt('Observations_Coreshift.txt').transpose()[0])
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
        MBH, M_UNIT, IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, CAM_SIZE_Y, FREQS_PER_DEC, FREQ_MIN, STEPSIZE = params_modelin
        f = open('model.in','r')
        text = f.readlines()
        f.close()
        text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
        text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
        text[6] = 'IMG_WIDTH\t(pixels)\t%d\n'%(IMG_WIDTH)
        text[7] = 'IMG_HEIGHT\t(pixels)\t%d\n'%(IMG_HEIGHT)
        text[8] = 'CAM_SIZE_X\t(Rg)\t\t%.15e\n'%(CAM_SIZE_X)
        text[9] = 'CAM_SIZE_Y\t(Rg)\t\t%.15e\n'%(CAM_SIZE_Y)
        text[10] = 'FREQS_PER_DEC\t(-)\t\t%d\n'%(FREQS_PER_DEC)
        text[11] = 'FREQ_MIN\t(Hz)\t\t%.15e\n'%(FREQ_MIN)
        text[13] = 'STEPSIZE\t(-)\t\t%.15e\n'%(STEPSIZE)
        f = open('model.in','w')
        f.writelines(text)
        f.close()

def Initialize_raptor_BHAC3D_modelc(params_raptor_BHAC_modelc, a):
    	# Set parameters in 'raptor_BHAC3D_model.c'
        Rhigh, Rlow = params_raptor_BHAC_modelc
        f = open('raptor_BHAC3D_model.c','r')
        text = f.readlines()
        f.close()
	text[47] = '\t\tsprintf(filename,"/\/vol/astro4/bhc_theory/jordyd/PhD/library/sane_a%s/data%%04d.blk",TIME_INIT);\n'%(a)
        if a == '+15o16':
		type = 'double'
	else:
		type = 'float'
	text[124] = '\t\tfseek(fp,-NPRIM*N1*N2*N3*sizeof(%s), SEEK_END);\n'%(type)
	text[126] = '\t\t%s buffer3;\n'%(type)
	text[142] = '\t\t\t\t\t\tfread(&buffer3,sizeof(%s),1,fp);\n'%(type)
	text[411] = '\treal Rhigh=%.15e;\n'%(Rhigh)
        text[412] = '\treal Rlow=%.15e;\n'%(Rlow)
        f = open('raptor_BHAC3D_model.c','w')
        f.writelines(text)
        f.close()

def Initialize_parametersh(params_parametersh):
    	# Set parameters in 'parameters.h'
    	LOG, TETRAD, VTKFILE, IMGFILE, SPECFILE, source_dist = params_parametersh
        f = open('parameters.h','r')
        text = f.readlines()
        f.close()
        text[31] = '#define LOG %d\n'%(LOG)
        text[32] = '#define TETRAD %d\n'%(TETRAD)
        text[34] = '#define VTKFILE (%d)\n'%(VTKFILE)
        text[35] = '#define IMGFILE (%d)\n'%(IMGFILE)
        text[36] = '#define SPECFILE (%d)\n'%(SPECFILE)
        text[137] = '#define source_dist    (%.15e) // Distance to M81 (cm)\n'%(source_dist)
        num_indices = num_Radio + num_IR + num_Coreshift
        text[229] = '#define num_indices %d\n'%(num_indices)
        f = open('parameters.h','w')
        f.writelines(text)
        f.close()

def Set_parameters():
    	# Load parameters from 'MultiNestRAPTORinput.txt'
    	_, params = np.loadtxt('MultiNestRAPTORinput.txt', dtype=str, delimiter='=', unpack=True)
	a = params[26]
    	params = params[0:26].astype(float)
       	# Set global variables
        global MBH_var, M_UNIT_var, Rhigh_var, i_var, Radio, IR, Coreshift
        MBH_var, M_UNIT_var, Rhigh_var, i_var, Radio, IR, Coreshift = params[0:7].astype(int)
        global MBH, M_UNIT, Rhigh, i, data_number
        if not MBH_var:
            	MBH = params[7]
        if not M_UNIT_var:
            	M_UNIT = params[8]
        if not Rhigh_var:
            	Rhigh = params[16]
        if not i_var:
            	i = params[24]
        data_number = params[25]
    	# Initialize input files RAPTOR
	Initialize_imgrendererc()
	Initialize_sizeflux_v4c(params[9], params[10], params[11], params[23])
    	Initialize_modelin(params[7:16])
    	Initialize_raptor_BHAC3D_modelc(params[16:18],a)
    	Initialize_parametersh(params[18:24])

def Set_modelin(MBH, M_UNIT):
        # Set MBH and M_UNIT in 'model.in'
        f = open('model.in','r')
        text = f.readlines()
        f.close()
        text[0] = 'MBH\t\t(g)\t\t%.15e\n'%(MBH)
        text[1] = 'M_UNIT\t\t(g)\t\t%.15e\n'%(M_UNIT)
        f = open('model.in','w')
        f.writelines(text)
        f.close()

def Set_raptorBHAC3D_modelc(Rhigh):
        # Set Rhigh in 'raptor_BHAC3D_model.c'
        f = open('raptor_BHAC3D_model.c','r')
        text = f.readlines()
        f.close()
        text[411] = '\treal Rhigh=%.15e;\n'%(Rhigh)
        f = open('raptor_BHAC3D_model.c','w')
        f.writelines(text)
        f.close()

def RAPTOR(MBH, M_UNIT, Rhigh, i, data_number):
        # Set parameter values for RAPTOR:
        Set_modelin(MBH, M_UNIT)
        Set_raptorBHAC3D_modelc(Rhigh)
    	# Run RAPTOR
    	os.system('make clean')
        os.system('make bhac CPU=1')
    	os.system('set OMP_NUM_THREADS=24')
    	os.system('./RAPTOR model.in %d %f'%(data_number, i))
    	# Load flux values
    	Freq, Flux = np.loadtxt('output/spectrum_%d_%.2f.dat'%(data_number,i), dtype=str, unpack=True)
	Freq = Freq.astype(float)
	Flux = Flux.astype(float)
    	# Write output to 'Run-name_flux.txt'
    	f = open("%s_Flux.txt"%(sys.argv[1]),'a+')
    	f.write(' '.join(map(str,Flux))+'\n')
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
    	if i_var: # i has a uniform distribution between 5 and 90
        	cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] = 85*cube[MBH_var+M_UNIT_var+Rhigh_var+i_var-1] + 5
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
    	Flux = RAPTOR(MBH, M_UNIT, Rhigh, i, data_number)
    	# Calculate chisquared values for radio, IR and Coreshift
    	Chisquare_Radio, Chisquare_IR, Chisquare_Coreshift = [0,0,0]
    	if Radio:
		_, Measurements_Radio, Sigma_Radio = np.loadtxt('Observations_Radio.txt').transpose()
		Measurements_Radio = np.multiply(Measurements_Radio, 1000)
		Sigma_Radio = np.multiply(Sigma_Radio, 1000)
        	for freqnum in range(num_Radio):
#            		Chisquare_Radio += (Measurements_Radio[freqnum] - np.log10(Flux[freqnum]))**2  /(2*Sigma_Radio[freqnum]**2)
			Chisquare_Radio += (Measurements_Radio[freqnum] - Flux[freqnum])**2  /(2*Sigma_Radio[freqnum]**2)
#			Chisquare_Radio += (np.log10(Measurements_Radio[freqnum]) - np.log10(Flux[freqnum]))**2 /(2*(np.log10(Measurements_Radio[freqnum])-np.log10(Measurements_Radio[freqnum] - Sigma_Radio[freqnum]))**2)
    	if IR:
		_, Measurements_IR = np.loadtxt('Observations_IR.txt').transpose()
		Measurements_IR = np.multiply(Measurements_IR, 1000)
                Sigma_IR = np.multiply(Measurements_IR, 0.2)
                for freqnum in range(num_IR):
                        Chisquare_IR += (Measurements_IR[freqnum] - Flux[num_Radio+freqnum])**2  /(2*Sigma_IR[freqnum]**2)
#			Chisquare_IR += (np.log10(Measurements_IR[freqnum]) - np.log10(Flux[num_Radio+freqnum]))**2 /(2*(np.log10(Measurements_IR[freqnum])-np.log10(Measurements_IR[freqnum] - Sigma_IR[freqnum]))**2)
    	if Coreshift:
#		_, Measurements_Major, Sigma_Major, Measurements_Minor, Sigma_Minor = np.loadtxt('Observations_Coreshift.txt').transpose()
		_, Measurements_Major, Err_up_Major, Err_down_Major = np.loadtxt('Observations_Coreshift.txt').transpose()
        	os.system('./sizeflux_v4_5 1 %d %.2f %f %d'%(num_Coreshift, i, MBH, data_number))
            	Wavelength, Major, Minor, _, _, _, _ = np.loadtxt('output/lambda_th_%d_%.2f.dat'%(data_number,i), dtype=float, delimiter=' ', unpack=True)
		for Wavelengthnum in range(len(Wavelength)):
                	Chisquare_Coreshift += (Measurements_Major[Wavelengthnum] - Major[Wavelengthnum])**2  /(2*Err_up_Major[Wavelengthnum]**2)
#			Chisquare_Coreshift += (Measurements_Minor[Wavelengthnum] - Minor[Wavelengthnum])**2  /(2*Sigma_Minor[Wavelengthnum]**2)
#			Chisquare_Coreshift += (np.log10(Measurements_Major[Wavelengthnum]) - np.log10(Major[Wavelengthnum]))**2 /(2*(np.log10(Measurements_Coreshift[Wavelenghtnum])-np.log10(Measurements_Coreshift[Wavelengthnum] - Err_down_Major[Wavelengthnum]))**2)
		# Write output to 'Run-name_Major.txt'
        	f = open("%s_Major.txt"%(sys.argv[1]),'a+')
        	f.write(' '.join(map(str,Major))+'\n')
        	f.close()
		# Write output to 'Run-name_Minor.txt'
		f = open("%s_Minor.txt"%(sys.argv[1]),'a+')
                f.write(' '.join(map(str,Minor))+'\n')
                f.close()
        # Calculate total chisquared value
    	Chisquare_Total = Chisquare_Radio/num_Radio + Chisquare_IR/num_IR + Chisquare_Coreshift/num_Coreshift
    	# Write Chisquared values to 'Run-name_Chisquared.txt'
    	f = open("%s_Chisquared.txt"%(sys.argv[1]),'a+')
    	f.write("%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n"%(MBH, M_UNIT, Rhigh, i, Chisquare_Radio/num_Radio,Chisquare_IR/num_IR,Chisquare_Coreshift/num_Coreshift,Chisquare_Total))
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
	os.system('python Convergenceplotter.py %s'%(sys.argv[1]))
    	# Return loglikelihood
    	return - Chisquare_Total

## Main code
Set_parameters()

parameters = np.array(["MBH", "M_UNIT", "Rhigh", "i"])[np.array([MBH_var,M_UNIT_var,Rhigh_var,i_var]).astype(bool)]
print(parameters)
n_params = len(parameters)
prefix = sys.argv[1]

result = pymultinest.run(LogLikelihood=myloglike, Prior=myprior, n_dims=n_params, importance_nested_sampling=False, multimodal=True, evidence_tolerance=0.5, sampling_efficiency=0.8, null_log_evidence=-1e90, max_modes=100, mode_tolerance=-1e90, outputfiles_basename=prefix, resume=True, verbose=True, write_output=True, n_live_points=400, max_iter=0)


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

