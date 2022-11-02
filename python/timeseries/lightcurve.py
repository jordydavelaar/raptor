import numpy as np
import matplotlib.pyplot as plt
from sys import argv

offset=0
avg = 0

N=int(argv[1])
start=int(argv[2])
dt = 10.

print(N)


frequencies=["" for x in range(1)]
units =["" for x in range(1)]

def plot_lightcurve(rh,inc):
        lightcurve=np.zeros((N,5))
        time=np.zeros(N)

        for t in range(start,start+N):
                file= "inc%d"%inc+"/R%d"%rh+"/output/spectrum_%d_%.02f.dat"%(t,inc)
                data = np.loadtxt(file)
                lightcurve[t-start]=data
                time[t-start]=(t-start)*dt
        lightcurve=np.transpose(lightcurve)
       	output_array = np.insert(lightcurve, 0, time, axis=0)

       	np.savetxt("total_lcurve_inc%d_R%d"%(inc,rh)+".dat",np.transpose(output_array))


rh=40
inc=50
plot_lightcurve(rh,inc)
