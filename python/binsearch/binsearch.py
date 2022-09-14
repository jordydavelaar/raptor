import numpy as np
from sys import argv

def BinSearch(index_min, index_max, index_step, M, inc, F_desire = 2.5):
    M = np.log10(M)

    if (len(M)==2):
        M_new = (M[0] + M[1])/2.
    
    else:
        F_list = []
        for index in range(index_min, index_max, index_step):
            f, F = np.loadtxt("output/spectrum_%d_%.02f.dat"%(index,inc), unpack=True, usecols = (0,1))
            
            if (np.isscalar(F) == True) and (f == 2.3e+11):
                F_list.append(F)

            else:
                F_list.append(F[np.where(f == 2.30e+11)[0][0]])
    
        #determine the average F over the dump files
        F_dummy = sum(F_list)/len(F_list)
        
        #and check after if everthing went according to plan
        F_file = open("F.txt", "a")
        F_file.write("{}\n".format(F_dummy))
        F_file.close()

        #read in M parameter file
        sorted_M = sorted(M)
        M_last = M[len(M)-1]
    
        if (F_dummy <= F_desire+0.01) and (F_dummy >= F_desire-0.01):
            M_new = M_last

        #we need to go left
        elif (F_dummy > F_desire):
            M_new = (M_last + sorted_M[np.where(sorted_M==M_last)[0][0]-1])/2.

        #we need to go right
        elif (F_dummy < F_desire):
            M_new = (M_last + sorted_M[np.where(sorted_M==M_last)[0][0]+1])/2.

    return 10**M_new

M_start = np.loadtxt("M.txt")
dump_min = int(argv[1])  #800
dump_max = int(argv[2])  #1000
dump_step = int(argv[3]) #10
inc = int(argv[4])

M_UNIT = BinSearch(dump_min, dump_max, dump_step, M_start,inc)

#write it to the txt, so we can keep track of what is used
file = open("M.txt", "a")
file.write("{}\n".format(M_UNIT))
file.close()

#print the current M_unit so one can give it as a parameter to raptor_multi.sh
print(M_UNIT)
