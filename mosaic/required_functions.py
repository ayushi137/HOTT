from __future__ import print_function
##############################################################################
'''
Functions need for preforming various modification on astronomical images.
Following are key functions:
- regrid (this uses cartesian)
- resconvolve
- reproject (This one doesn't always works)

For more information, read the discription for each function. 

By: Ayushi Singh
Last Modified: November 16, 2019
'''
##############################################################################

# package required 
from astropy.io import fits
from astropy import wcs
import numpy as np 
from scipy import signal
import scipy.interpolate
import time, datetime


def regrid(sourceimage,targetimage, savefile = True, resultimage = 'same', header = None, mode = 3):
	"""
    Regrid sourceimage to the size and shape of targetimage. 
    This uses the another function called cartesian. Read in cartesian for more information. 

    Parameters
    ----------
    sourceimage : string 
    	path to image to regrid - should already be convolved to appropriate 
        resolution using resconvolve.
    targetimage : string
    	path to image whose grid is to be used in regridding.
    savefile : boolean
    	If True then saves the file using the parameter resultimage and header.
    	If False then returns the regird image.
    	Default is True.	
    resultimage : string
    	If 'same', then use the same directoy and file name is sourceimage.resamp.fits.
    	If specified then use the given path and file name.
    	Default is 'same'.
    header : header format or None
    	If None, then The header of targetimage is used. But the wavelength of sourceimage is used (if given). 
    	If heater is given, then that is used for result header.
    	Default is None. 
    mode : integer 
    	Degrees of the bivariate spline.
    	Default is 3, which is cubic.	

    Returns
    -------
    if savefile is True
    saves a file : fits
        stores an 2-D image of sourceimage regrid to targetimage at resultimage. 
    if savefile is False
    tofill : 2-D array 
    	sourceimage regird to the shape and size of targetimage	      

    ~~~~~ By: Natalie Price-Jones ~~~~~
   	"""
	# Start timer
	start = time.time()
	# Load in source data and header information
	sdata,sheader = fits.getdata(sourceimage,header = True)
	# nan to zero
	sdata[np.isnan(sdata)] = 0
	sdata[np.isinf(sdata)] = 0
	# Create WCS object for the source image
	sourcewcs = wcs.WCS(sheader)
	# Create array of pixel indices in source image
	x = np.arange(sdata.shape[1])
	y = np.arange(sdata.shape[0])
	# Interpolate the image data over pixel indices
	interp = scipy.interpolate.RectBivariateSpline(y,x,sdata, kx=mode, ky=mode)
	# Load in target grid data
	tdata,theader = fits.getdata(targetimage,header=True)
	tdata[np.isnan(tdata)] = 0
	# Create WCS object for target grid
	targetwcs = wcs.WCS(theader)
	# Create all possible pairs of pixel coordinates in target grid
	coords = cartesian([np.arange(tdata.shape[1]),np.arange(tdata.shape[0])])
	# Extract x and y columns of pixel pairs
	xpixs = coords[:,0]
	ypixs= coords[:,1]
	# Convert target grid pixels to ra/dec 
	world = targetwcs.wcs_pix2world(coords,0)
	# Convert target grid ra/dec to source pixel coordinates
	dpix = sourcewcs.wcs_world2pix(world,0)
	# Extract x and y columns of converted pixel pairs
	xdpixs = dpix[:,0]
	ydpixs = dpix[:,1]
	# Find where target grid corresponds to actual source image data
	good = np.where((xdpixs >= min(x)) & (xdpixs <= max(x)) & (ydpixs >= min(y)) & (ydpixs <= max(y)))
	# Pick out only pixels with relevant image data
	xpixs = xpixs[good]
	ypixs = ypixs[good]
	xdpixs = xdpixs[good]
	ydpixs = ydpixs[good]
	# Create grid to fill up with source image regrid
	tofill = np.copy(tdata)
	tofill[:] = np.NAN
	# Loop over relevant pixels and fill up source image regrid
	for i in range(len(ypixs)):
	    ypix = ypixs[i]
	    xpix = xpixs[i]
	    xdpix = xdpixs[i]
	    ydpix = ydpixs[i]
	    tofill[ypix,xpix] = interp(ydpix,xdpix)
	# End timer
	end = time.time()
	# Print time to run
	print ('Regridded in ', str(datetime.timedelta(seconds=end-start)))
	if savefile == True:
		# name of the path based on the resultimage
		if resultimage == 'same':
			path = sourceimage[:-5]+'.resamp.fits'
		else:
			path = resultimage	

		# change the wavelength if it's there
		if header == None:
			try:
				theader['WAVELNTH'] = (sheader['WAVELNTH'], '[micrometer] The reference wavelength' ) 
			except:
				print ("Wavelength is not given in the header")
		else:
			theader = header			      

		# saves as the fits file
		print ("File Saved as:", path)	
		fits.writeto(path, tofill, theader)

		print ("Regriding done. File is saved. \n")

		return
		
	elif savefile == False:
		print ("Regriding done. Image is returned. \n")
		return tofill	 


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.
    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.
    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

	from: https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/extmath.py
    """
    arrays = [np.asarray(x) for x in arrays]
    shape = (len(x) for x in arrays)
    dtype = arrays[0].dtype

    ix = np.indices(shape)
    ix = ix.reshape(len(arrays), -1).T

    if out is None:
        out = np.empty_like(ix, dtype=dtype)

    for n, arr in enumerate(arrays):
        out[:, n] = arrays[n][ix[:, n]]

    return out

def resconvolve(sourceimage, highres, lowres, resultimage = 'same'):
	"""
    Convolve sourceimage to a lower resolution image. 

    Parameters
    ----------
    sourceimage : string 
    	Rath to image to be convolved.
    highres : float
    	Resolution of the current image (in arcsec).
    lowres : float
    	Resolution of the requied convolved image (in arcsec).
    resultimage : string
    	If 'same', then use the same directoy and file name is sourceimage.resamp.fits.
    	If specified then use the given path and file name.	
    	Default is 'same'.

    Returns
    -------
    saves a file : fits
        Stores an 2-D image of sourceimage convolved to the required resolution at resultimage. 
        The header of sourceimage is used.

   	"""
   	# Start timer
	start = time.time()
	# Load in source data and header information
	data,header = fits.getdata(sourceimage,header = True)
	# nan to zero
	data[np.isnan(data)] = 0
	data[np.isinf(data)] = 0
	# get the pixelsize in arcsec, where the value is in degree*3600(arcsec/degree)
	pixelsize = header['CDELT2']*3600
	# gaussian consant to convert sigma to FWHM 2*sqrt(2ln(2))
	constant = 2*np.sqrt(2*np.log(2))
	# FWHM of current image 
	FWHM_low = lowres/pixelsize
	# sigma of current image
	sigma_low = FWHM_low/constant
	# FWHM of required resolution
	FWHM_high = highres/pixelsize
	# sigma of required resolution
	sigma_high = FWHM_high/constant
	# FWHM calulated for the gaussian of convolution kernel
	FWHM = np.sqrt(FWHM_low**2 - FWHM_high**2)
	# sigma for the gaussian of convolution kernel
	sigma = FWHM/constant
	# number of x and y pixels
	x = int(FWHM)*10+1
	y = int(FWHM)*10+1
	# making the 2-D image of the convolution kernel by making 2, 1-D gaussian
	# and normalized by the gaussian normalization factor
	gauss1 = signal.general_gaussian(x,1, sigma)/((sigma)*np.sqrt(2*np.pi))
	gauss2 = signal.general_gaussian(y,1, sigma)/((sigma)*np.sqrt(2*np.pi))
	# combine these two to create 2D image
	kernel = np.outer(gauss1, gauss2)
	# convolve the image using signal.fftconvolve premade function and kernel
	convolved = signal.fftconvolve(data, kernel, mode='same')

	# change the resolution if it's there
	header['RESO'] = (lowres, '[arcsec] Resolution of the map')

	# name of the path based on the resultimage
	if resultimage == 'same':
		path = sourceimage[:-5]+'.conv.fits'
	else:
		path = resultimage

	# End timer
	end = time.time()
	# Print time to run
	print ('Convolving in ', str(datetime.timedelta(seconds=end-start)))
	
	# saves as the fits file
	print ("File Saved as:", path)		
	fits.writeto(path, convolved, header)

	
	print ("Convolving done. File is saved. \n")


	return 

def reproject(sourceimage,targetimage, savefile = True, resultimage = 'same', header = None, orders = 3):
	"""
    Regrid sourceimage to the size and shape of targetimage. 
    This uses the another function called cartesian. Read in cartesian for more information. 

    Parameters
    ----------
    sourceimage : string 
    	path to image to regrid - should already be convolved to appropriate 
        resolution using resconvolve.
    targetimage : string
    	path to image whose grid is to be used in regridding.
    savefile : boolean
    	If True then saves the file using the parameter resultimage and header.
    	If False then returns the regird image.
    	Default is True.	
    resultimage : string
    	If 'same', then use the same directoy and file name is sourceimage.resamp.fits.
    	If specified then use the given path and file name.
    	Default is 'same'.
    header : header format or None
    	If None, then The header of targetimage is used. But the wavelength of sourceimage is used (if given). 
    	If heater is given, then that is used for result header.
    	Default is None. 
    order : integer 
		Order of interpolation.	
		Default 3	

    Returns
    -------
    if savefile is True
    saves a file : fits
        stores an 2-D image of sourceimage regrid to targetimage at resultimage. 
    if savefile is False
    newdata : 2-D array 
    	sourceimage regird to the shape and size of targetimage	      

   	"""
	from reproject import reproject_interp
	# Start timer
	start = time.time()
	# Load in source data and header information
	sdata = fits.open(sourceimage)[0]
	sheader = sdata.header
	tdata = fits.open(targetimage)[0]
	theader = tdata.header
	newdata, footprint = reproject_interp(sdata, theader, order=orders) # 3 = bicubic
	# End timer
	end = time.time()
	# Print time to run
	print ('Regridded in ', str(datetime.timedelta(seconds=end-start)))

	# change the wavelength if it's there
	if header == None:
		try:
			theader['WAVELNTH'] = (sheader['WAVELNTH'], '[micrometer] The reference wavelength' ) 
		except:
			print ("Wavelength is not given in the header")
	else:
		theader = header			
	
	if savefile == True:
		# name of the path based on the resultimage
		if resultimage == 'same':
			path = sourceimage[:-5]+'.resamp.fits'
		else:
			path = resultimage			      

		# saves as the fits file	
		print ("File Saved as:", path)	
		fits.writeto(path, newdata, theader)

		print ("Regriding done. File is saved. \n")

		return

	elif savefile == False:
		print ("Regriding done. Image is returned. \n")
		return newdata, theader

