import numpy as np
import h5py
import click
import sys,os,glob
import itertools as it

programDir = os.path.dirname(os.path.realpath(__file__))

@click.command()
@click.option("--dumpfile",default = "/dev/null")
@click.option("--outputdir",default = programDir)
def convertFileForRaptor(dumpfile,outputdir):
    hfp = h5py.File(dumpfile,'r')
    with open(os.path.join(outputdir,dumpfile.split('/')[-1].replace(".h5",".txt")),"w") as fp:
        header = hfp['header']
        N1 = header['n1'][()]
        N2 = header['n2'][()]
        N3 = header['n3'][()]
        gam = header['gam'][()]
        a = header['a'][()]
        dx1,dx2,dx3 = header['dx1'][()],header['dx2'][()],header['dx3'][()]
        startx1,startx2,startx3 = header['startx1'][()],header['startx2'][()],header['startx3'][()]
        hslope = header['hslope'][()]
        Rin = header['r_in'][()]
        Rout = header['r_out'][()]
        mks_smooth = header['mks_smooth'][()]
        poly_xt = header['poly_xt'][()]
        poly_alpha = header['poly_alpha'][()]


        fp.write(f"{N1} {N2} {N3} {gam} {a} {dx1} {dx2} {dx3} {startx1} \
            {startx2} {startx3} {hslope} {mks_smooth} {poly_xt} \
                {poly_alpha} {Rin} {Rout} \n")
        prims = hfp['prims']
        nprims = hfp['header']['np'][()]
        for i,j,k in it.product(range(N1),range(N2),range(N3)):
            for l in range(nprims):
                fp.write(f"{prims[i,j,k,l]} ")
            fp.write("\n")

if __name__ == "__main__":
    convertFileForRaptor()
