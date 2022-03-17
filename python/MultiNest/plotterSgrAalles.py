import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys



# argv
incl = float(sys.argv[1])
prefix = sys.argv[2]
#MBH = sys.argv[3]
#M_UNIT = sys.argv[4]
#a = sys.argv[5]
datanumber = 900

# Load data
os.chdir(".")
Radio_freq, Radio_flux, Radio_err = np.loadtxt('Observations_Radio.txt').transpose()
Radio_flux = np.multiply(Radio_flux, 1000)
IR_freq, IR_flux = np.loadtxt('Observations_IR.txt').transpose()
IR_flux = np.multiply(IR_flux, 1000)
Size_wavelength, Size_size, Size_err_up, Size_err_down = np.loadtxt('Observations_Coreshift.txt').transpose()


# Plot spectrum
fig = plt.figure(num=None, figsize=(10,4))
plt.subplot(1,2,1)

plt.errorbar(Radio_freq, Radio_flux, yerr=Radio_err, fmt='.',color='k', label='Radio')
plt.scatter(IR_freq, IR_flux, label='IR', color='k', marker='v')

# Load spectrum simulation
Simulation_freq, Simulation_flux = np.loadtxt('output/spectrum_%d_%.2f.dat'%(datanumber,incl)).transpose()

plt.plot(Simulation_freq, Simulation_flux, color='b')

xlabelsize = '20'
ylabelsize = '20'
titlesize = '25'

plt.xlabel('$\\nu$ (Hz)', size=xlabelsize)
plt.ylabel('$F_{\\nu}$ (mJy)', size=ylabelsize)
plt.xlim(1e10, 1e16)
plt.ylim(1e-2, 1e5)
plt.xscale('log')
plt.yscale('log')
#plt.title('Spectrum M81* frame=%d i=%.2f'%(frame,incl), size=titlesize)
#plt.title('a = %s: MBH = %s g; M_UNIT = %s g; i = %s'%(a, MBH, M_UNIT, incl), size=15)
#plt.legend(loc=0)
#plt.savefig('Spectrum/Spectrum_%d_%.2f.png'%(frame,incl))

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)


# Core size plot
plt.subplot(1,2,2)


Size_freq = np.divide(29979245800, Size_wavelength)

plt.errorbar(Size_freq, Size_size, yerr=Size_err_up, fmt='.', color='k',label='Observations') # eigenlijk wil je error omhoog en omlaag apart 


Simulation_wavelength, Simulation_size, _, _, _, _, _ = np.loadtxt('output/lambda_th_%d_%.2f.dat'%(datanumber,incl)).transpose()

Simulation_freq = np.divide(29979245800, Simulation_wavelength)


plt.plot(Simulation_freq, Simulation_size, color='b')

xlabelsize = '20'
ylabelsize = '20'
titlesize = '25'

plt.xlabel('$\\nu$ (Hz)', size=xlabelsize)
plt.ylabel('Core size (mas)', size=ylabelsize)
plt.xlim(1e10, 1e12)
plt.ylim(1e-2, 1e1)
plt.xscale('log')
plt.yscale('log')
#plt.title('Spectrum M81* frame=%d i=%.2f'%(frame,incl), size=titlesize)
#plt.title('a = %s: MBH = %s g; M_UNIT = %s g; i = %s'%(a, MBH, M_UNIT, incl), size=15)
#plt.legend(loc=0)
#plt.savefig('Spectrum/Spectrum_%d_%.2f.png'%(frame,incl))

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)





fig.tight_layout()
#fig.savefig('M81data/Spectrum/Spectrum_%d_%.2f.png'%(frame,incl))
#fig.savefig('M81data/Spectrum/Bestfitspectrum%d.png'%(frame))
fig.savefig('Spectrum_%s.png'%(prefix))
