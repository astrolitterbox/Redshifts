#rewrite of a Perl script by D. Mast and S. Sanchez
import numpy
import pyfits

def z_to_factor(z_in, z_out):
	#calculate ratio of linear dimensions

#def luminosity_distance(z_in, z_out):



def flux_scaling(cube, z_in, z_out):
	#multiply the values of binned spaxels by the flux ratio
	2.5*log(f_o/f_i) = (luminosity_distance_in/

def redshift_distance(z_in, z_out):
	rel_z = (1 + z_out)/(1 + z_in) - 1
	return rel_z

def binCube(cube, factor, outputFileName):
	#rebin the cube
	#multiply by ratio of sum of pixel values before/after binning


def cosmologicalDimming(cube, z_in, z_out):
	


def spectrum_shifting(cube, z_in, z_out):
	#slice the cube, lambda = the spectral dimension
	lambda = lambda*(1 + redshift_distance(z_in, z_out)
	#slice the cube again -- keep only the relevant wavelength range
	


#etc:
#different psfs
#different instrument sensitivities: sky background

def main():
  inputFileName = ''
  inputFile = pyfits.open(inputFileName)
  inputImage = inputFile[0].data

