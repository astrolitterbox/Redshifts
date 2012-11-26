#rewrite of a Perl script by D. Mast and S. Sanchez
import numpy as np
import pyfits
import math 
import distance as cosmology
import matplotlib.pyplot as plt
from califa_cmap import califa_vel, califa_int
from congrid import congrid
import make_2d_gaussian as make_gauss
import scipy.signal
# -------------------------- A few ideas:
#1: Surface brightness dimming
#2: Set exposure time as free variable, depending on S/N
#3: PSF
#4: Various binning scenarios

fwhm = 3.7
z_out = 0.05
z_in = 0.018
psf = make_gauss.make_2d_gaussian(17, fwhm)
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
	return cube/flux_ratio


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

def get_convolved_image(cube, psf):
	return scipy.signal.convolve2d(cube,psf,mode='same',boundary='fill', fillvalue=0)
	

def get_deconvolved_image(cube, psf):
	return scipy.signal.deconvolve(cube, psf)
		
#multiply by ratio of sum of pixel values before/after binning


#def cosmologicalDimming(cube, z_in, z_out):
	
#def spectrum_shifting(cube, z_in, z_out):
	#slice the cube, lambda = the spectral dimension
	#lambda = lambda*(1 + redshift_distance(z_in, z_out)
	
	#slice the cube again -- keep only the relevant wavelength range
	


#etc:
#different psfs
#different instrument sensitivities: sky background


filename = 'data/ARP220.COMB.rscube.fits '


fitsdata,fitshdr=read_fits_img(filename)


factor = get_magnification(z_in, z_out)
cube = fitsdata

print cube.shape
rebinned = np.empty((cube.shape[0], int(round(cube.shape[1]*factor, 0)), int(round(cube.shape[2] * factor, 0))))
print rebinned.shape, 'shape of resized cube', 1/factor**2, '/factor**2'

outputShape = int(round(cube.shape[1]*factor, 0)), int(round(cube.shape[2] * factor, 0))

#for i in range(0, cube.shape[0]):
for i in range(1450, 1451):
	cube_slice = cube[i, :, :]
	cube_slice = get_deconvolved_image(cube_slice, psf)
	rebinned[i, :, :] = congrid(cube_slice, outputShape, method='linear', centre=True, minusone=False)
	rebinned[i, :, :] = (1/factor**2)*rebinned[i, :, :]
	rebinned[i, :, :] = get_scaled_flux(rebinned[i, :, :], z_out)
	rebinned[i, :, :] = get_convolved_image(rebinned[i, :, :], psf)



print np.sum(rebinned), np.sum(cube)
#h=pyfits.PrimaryHDU(rebinned)
#hdu2=pyfits.HDUList([h])
#if os.path.exists('convolved.fits'): os.unlink('convolved.fits')
#hdu2.writeto('resized.fits')

plot_img(rebinned, 'rebinned')
#filename = 'resized.fits'
#fitsdata,fitshdr=read_fits_img(filename)
#plot_img(fitsdata[1450,:, :],'res_'+filename[:-4])


#plot_img((1/factor**2)*rebinned, 'reb'+filename[5:-4])

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

