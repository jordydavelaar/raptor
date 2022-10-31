
import numpy as np
from astropy.io import fits

from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time as aTime
import os
import sys


#constants
MSUN = 1.989e33
KPC  = 3.086e21

G  = 6.674e-8
c2 = 8.988e20

asindeg   = 206264.806
DEGPERMAS = 2.778e-7

#user should change these...
M  = 4e6* MSUN
d  = 8 * KPC

a    = 0.9375
R    = 40
inc  = 30
freq = 230e9

fov    = 40. #Field of view in Rg
pixels = 200 #pixels in file

#computing some necessary quantities
mas = ((G*M/c2)/d )* asindeg * 1000
dr = fov/pixels
PIX_SIZE = dr * mas

#passed to command line
ind = int(sys.argv[1])

def read_image(folder,ind,freq,inc):
    data =np.loadtxt("output/uniform_img_%.02e_%d.dat"%(freq,ind),skiprows=0,usecols=[2,3,4,5,6,7],unpack=True)
    print(data.shape)
    image=np.reshape(data,(6,pixels,pixels))
    image= np.transpose(image,axes=[0,2,1])
#    image = image[:,::-1,:]
    return image

def write_fits(folder,image,ind, a, R, inc, freq):

    user = "BHAC group"
    source = "Sgr A*"
    date = "2022-08-10"
    metric = "MKS"

    history = "metric=%s a=%f Rhigh=%f"%(metric, a, R)
    #get source postion
    loc=SkyCoord.from_name(source)

    #get RA and DEC in degree
    ra=loc.ra.deg
    dec=loc.dec.deg

    #convert date to mjd (modified julian date)
    modjuldate=(aTime(date)).mjd

    outfile = '.fits'

    header = fits.Header()
    header['AUTHOR']   =  user
    header['OBJECT']   =  source
    header['CTYPE1']   = 'RA---SIN'
    header['CTYPE2']   = 'DEC--SIN'
    header['CDELT1']   =  -PIX_SIZE*DEGPERMAS
    header['CDELT2']   =  PIX_SIZE*DEGPERMAS
    header['OBSRA']    = ra
    header['OBSDEC']   = dec
    header['FREQ']     = freq
    header['MJD']      = float(modjuldate)
    header['TELESCOP'] = 'VLBI'
    header['BUNIT']    = 'JY/PIXEL'
    header['STOKES']   = 'IQUVtautauf'
    header['HISTORY']  = history
    hdu = fits.PrimaryHDU(image,header=header)
    hdulist = [hdu]

    hdulist = fits.HDUList(hdulist)

    # Check whether the specified path exists or not
    if not os.path.exists(folder):
        os.makedirs(folder)

    hdu.writeto(folder+"/data%d.fits"%ind,overwrite=True)

image=read_image("",ind,freq,inc)
write_fits("fits/",image,ind,a,R,inc,freq)
