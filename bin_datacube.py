#rewrite of a Perl script by D. Mast and S. Sanchez
import numpy as np
import pyfits
import math 
import distance as cosmology
import matplotlib.pyplot as plt
from califa_cmap import califa_vel, califa_int
from congrid import congrid
# -------------------------- A few ideas:
#1: Surface brightness dimming
#2: Set exposure time as free variable, depending on S/N
#3: PSF
#4: Various binning scenarios



cosmo = {'omega_M_0':0.272, 'omega_lambda_0':0.728, 'omega_k_0':0.0, 'h':0.7}

cosmo = cosmology.set_omega_k_0(cosmo)

def get_magnification(z_in, z_out):
	#calculate ratio of linear dimensions
	lumdist_in = cosmology.luminosity_distance(z_in, **cosmo)	
	lumdist_out = cosmology.luminosity_distance(z_out, **cosmo)
#	print lumdist_in
	magnification = (lumdist_in/lumdist_out) *((1.+z_out)**2/(1.+z_in)**2)         #*(p_hi/p_lo))[0] -- we assume pixel scale is the same
	return magnification

def get_scaled_flux(cube, z_out):
	flux_ratio = (1+z_out)**4
	print flux_ratio, 'flux ratio'
	return cube/flux_ratio
	
def rebin(data, new_shape):
    N = 4
    height, width = data.shape
    #data = np.average(np.split(np.average(np.split(data, width // N, axis=1), axis=-1), height // N, axis=1), axis=-1)
    

### START rfits_img
def read_fits_img(filename):

    # READ FITS FILE
    fitsdata=pyfits.getdata(filename)
    fitshdr=pyfits.getheader(filename)
    nx = fitshdr['NAXIS1']
    ny = fitshdr['NAXIS2']
    #out=np.zeros((ny,nx))
    #infinite=np.isfinite(fitsdata,out)
    #fitsdata=fitsdata*out
    fitsdata=np.nan_to_num(fitsdata)
    return fitsdata,fitshdr
    
    
def plot_img(fitsdata,filename):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    cax = ax1.imshow(fitsdata, interpolation='nearest', cmap=califa_vel) 
    ax1.set_title(filename)
    cbar = fig1.colorbar(cax)
    plt.savefig(filename)

#def get_convolved_image(cube, 


		
#multiply by ratio of sum of pixel values before/after binning


#def cosmologicalDimming(cube, z_in, z_out):
	
#def spectrum_shifting(cube, z_in, z_out):
	#slice the cube, lambda = the spectral dimension
	#lambda = lambda*(1 + redshift_distance(z_in, z_out)
	
	#slice the cube again -- keep only the relevant wavelength range
	


#etc:
#different psfs
#different instrument sensitivities: sky background


filename='data/ARP220.V1200.rscube.fits'


fitsdata,fitshdr=read_fits_img(filename)

#print fitsdata[1450,:, :].shape

plot_img(fitsdata[1450,:, :],filename[5:-4])

factor = get_magnification(0.015, 0.05)
#print factor, 'factor'
cube = fitsdata[1450,:, :]

#factor = 0.5
#print int(round(cube.shape[0]*factor, 0)), int(round(cube.shape[1] * factor, 0)), 'new shape'
#cube = get_scaled_flux(cube, 0.05)
#cube = np.ones((4, 4))
#cube[1:3, 1:3] = 5
#cube[2, 2] = 8


rebinned = congrid(cube, (int(round(cube.shape[0]*factor, 0)), int(round(cube.shape[1] * factor, 0))), method='linear', centre=True, minusone=False)

 #(cube.shape[0]/rebinned.shape[0])*(cube.shape[1]/rebinned.shape[1])

print np.sum((1/factor**2)*rebinned), np.sum(cube)
plot_img((1/factor**2)*rebinned, 'reb'+filename[5:-4])

exit()

redshifts = 0.01*np.arange(1, 5)
size = np.empty((redshifts.shape))
for i, x in enumerate(redshifts):
	print i
	size[i] = z_to_factor(0.002, x)
print size



fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.set_yscale('log')

p = plt.scatter(redshifts, size, s=10, edgecolor="None", marker='o')
plt.savefig('sizes')

