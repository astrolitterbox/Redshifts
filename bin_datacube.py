#rewrite of a Perl script by D. Mast and S. Sanchez
import numpy as np
import pyfits
import math 
import distance as cosmology
import matplotlib.pyplot as plt
from califa_cmap import califa_vel, califa_int


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

def flux_scaling(cube, z_in, z_out):
	print 'kaciukai'
	#multiply the values of binned spaxels by the flux ratio
	#2.5*log(f_o/f_i) = (luminosity_distance_in/luminosity_distance_in)^2



def binCube(cube, factor):
	rebinnedCube = np.resize(cube, factor*cube.shape)
	return rebinnedCube

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




		
#multiply by ratio of sum of pixel values before/after binning


#def cosmologicalDimming(cube, z_in, z_out):
	
#def spectrum_shifting(cube, z_in, z_out):
	#slice the cube, lambda = the spectral dimension
	#lambda = lambda*(1 + redshift_distance(z_in, z_out)
	
	#slice the cube again -- keep only the relevant wavelength range
	


#etc:
#different psfs
#different instrument sensitivities: sky background


print get_magnification(0.002, 0.05)
print get_magnification(0.00331, 0.05)
filename='data/ARP220.V1200.rscube.fits'


fitsdata,fitshdr=read_fits_img(filename)

print fitsdata[1450,:, :].shape

plot_img(fitsdata[1450,:, :],filename[5:-4])

factor = get_magnification(0.015, 0.05)
rebinned = binCube(fitsdata[1450,:, :], factor)
plot_img(rebinned, 'reb'+filename[5:-4])

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

