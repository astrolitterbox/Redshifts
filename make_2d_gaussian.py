#
# make up a 2d gaussian
#
import numpy, pyfits, os, math

def gaussian_ftn(y,x):
    global xc, yc, sigma
    xterm = ((x-xc)/sigma)**2 
    yterm = ((y-yc)/sigma)**2 
    value =  math.exp(-0.5*(xterm+yterm))
    return value

def make_2d_gaussian(size, fwhm):
	global xc, yc, sigma
	xc = float(size-1)/2.
	yc = float(size-1)/2.
	sigma = fwhm/2.35482

	psf = numpy.zeros((size,size), dtype='f')
	for ix in range(size):
	  for iy in range(size):
	    x1 = float(ix) - 0.5
	    y1 = float(iy) - 0.5
	# crude integration by taking a 3x3 average within the pixel
	    val2 = ( gaussian_ftn(y1+0.1667,x1+0.1667) + gaussian_ftn(y1+0.1667,x1+0.5) +  gaussian_ftn(y1+0.1667,x1+0.8333) + \
		     gaussian_ftn(y1+0.5,x1+0.1667) + gaussian_ftn(y1+0.5,x1+0.5) +  gaussian_ftn(y1+0.5,x1+0.8333) + \
		     gaussian_ftn(y1+0.8333,x1+0.1667) + gaussian_ftn(y1+0.8333,x1+0.5) +  gaussian_ftn(y1+0.8333,x1+0.8333) )
	    psf[iy][ix] = val2/9.0

	psf_sum = psf.sum() 
	psf = psf/psf_sum   # normalise total value to 1
	h=pyfits.PrimaryHDU(psf)
	#hdu2=pyfits.HDUList([h])
	#if os.path.exists(outputFilename+'.fits'): os.unlink(outputFilename+'.fits')
	h.writeto('output/psf_new.fits', clobber=True)
	return psf


#make_2d_gaussian(51, 3.7, 'psf')	
