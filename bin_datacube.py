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
import sys

# -------------------------- A few ideas:
#1: Signal to noise degradation: should I just add Poisson noise after convolution?
#2: PSF: deconvolution
#3: redshifting: changes in spectral sampling (resolution)
#4: Various binning scenarios
#2: Set exposure time as free variable, depending on S/N

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}

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
	

def get_deconvolved_image(image, psf):
	#OriginalImage = inverse_FFT[ FFT(ObservedImage) / FFT(PointSpreadFunction) ]
	deconvolved_image = np.fft.irfft2(scipy.signal.wiener((np.fft.rfft2(image))/scipy.signal.wiener(np.fft.rfft2(psf))))
	#h=pyfits.PrimaryHDU(deconvolved_image)
	#hdu2=pyfits.HDUList(h)
	#hdu2.writeto('deconv_'+'.fits')
	return deconvolved_image
		


def do_redshifting(cube, factor, psf, convolve=True):
  outputShape = int(round(cube.shape[1]/factor, 0)), int(round(cube.shape[2]/factor, 0))
  rebinned = np.empty((cube.shape[0], int(round(cube.shape[1]/factor, 0)), int(round(cube.shape[2] / factor, 0))))
  print cube.shape, 'cshape', psf.shape, 'psf', rebinned.shape, 'rebinned', outputShape, 'os'
  #fluxFactor = (cube.shape[1]/rebinned.shape[1])*(cube.shape[2]/rebinned.shape[2])

  #print outputShape, 'os',  areaFactor, 'area scale'

  for i in range(0, cube.shape[0]):
	  cube_slice = cube[i, :, :]
	  if convolve == True:
	    cube_slice = get_convolved_image(cube_slice, psf)
	  #rebinned = np.real(get_deconvolved_image(cube_slice, psf))
	  #rebinned[i, :, :] = get_scaled_flux(rebinned[i, :, :], z_out)	  
	  rebinned[i, :, :] = congrid(cube_slice, outputShape, method='linear', centre=True, minusone=False)
	  areaFactor=(cube.shape[1]*cube.shape[2])/(rebinned.shape[1]*rebinned.shape[2])
	  rebinned[i, :, :] = areaFactor*rebinned[i, :, :]


  return rebinned



#etc:
#different psfs
#different instrument sensitivities: sky background


galaxy_name = sys.argv[1]

filename = 'data/Z0_'+galaxy_name+'.rscube.fits'



z_out = 0.05 
z_in = 0.01728 #0.002392  for NGC1637 #0.001728 #for NGC1058, 



# ----------------------- data extension -- 
hdulist = pyfits.open(filename)
#hdulist.info()
fitsdata = hdulist[0].data
fitshdr = hdulist[0].header
cube = fitsdata
#factor = get_magnification(z_in, z_out)
factor = 3
fwhm = 3.5*factor
psf = make_gauss.make_2d_gaussian(101, fwhm)
#psf = psf[3:-3, 1:]
print psf.shape, 'psf'
print factor, 'factor'
rebinned = do_redshifting(cube, factor, psf, True)
print rebinned.shape, 'shape of resized cube', 1/factor**2, '/factor**2'
print np.sum(rebinned), np.sum(cube)
#rebinned = np.swapaxes(rebinned, 1, 2)
h=pyfits.PrimaryHDU(rebinned, header=fitshdr)
#h=pyfits.PrimaryHDU(cube, header=fitshdr)
hdu2=pyfits.HDUList([h])
#if os.path.exists('convolved.fits'): os.unlink('convolved.fits')
hdu2.writeto('output/Z1_'+galaxy_name+'_conv_before.fits', clobber=True)

exit()

# ----------------------- errors ------------------


hdulist = pyfits.open(filename)
hdulist.info()
fitsdata = hdulist['ERROR'].data
fitshdr = hdulist[0].header
cube = fitsdata

factor = get_magnification(z_in, z_out)
print factor, factor*(cube.shape[1]), factor*(cube.shape[2])
print cube.shape, 'cube shape'
#rebinned = do_redshifting(cube, factor, psf, False)
#rebinned = rebinned*factor**2
#print rebinned.shape, 'shape of resized cube', 1/factor**2, '/factor**2'
#print np.sum(rebinned), np.sum(cube)
#h=pyfits.PrimaryHDU(rebinned, header=fitshdr)
h=pyfits.PrimaryHDU(cube, header=fitshdr)
hdu2=pyfits.HDUList([h])
#if os.path.exists('convolved.fits'): os.unlink('convolved.fits')
hdu2.writeto('NGC4676B_err_big.fits', clobber=True)












exit()
#plot_img(cube[1450,:, :], filename[5:-21]+'_deconv_slice')
#filename = 'resized.fits'

#plot_img(fitsdata[1450,:, :], filename[5:-21])


#plot_img(rebinned[1450,:,:], 'conv_deconv_'+filename[5:-21])
#plot_img(rebinned[1450, :, :], 'conv_'+filename[5:-21])
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

