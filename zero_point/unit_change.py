from __future__ import print_function
##############################################################################
'''
Reads a herschel fits file and change the unit to MJy/sr. Also save extract 
the sigma and coverage fit files and move them to appropriate folder.

This script requires following functions and files:
initial_parameter.py

To run:
run unit_change.py 'Field'

Example:
run unit_change.py TauS2

By: Ayushi Singh
Last Modified: 30 December 2019
'''
##############################################################################

from astropy.io import fits
import numpy as np
import os
import sys
import glob

# importing other scipts
from initial_parameter import path_to_folder, beamsize

def unit_change(file_name,field, wavelength, resultimage = 'same'):
	"""
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    file_name : string
        path to image. 
    field : string
    	name of the field    
    wavelength : float
    	wavelength of the image     
    resultimage : string
    	if 'same', then use the same directoy and file name is sourceimage.resamp.fits.
    	if specified then use the given path and file name.

    Returns
    -------
    saves 3 maps : fits
        the three files are data, sigma and coverage map
        stores an 2-D image of file_name with proper units in MJy/sr. 
        The header of file_name is used with correctly altered.  

    """
    # read the data and the header
	all_data = fits.open(file_name)
	data = all_data[1].data
	header = all_data[1].header
	header_true = all_data[0].header
	

	# multiply the data by a constant to change the units
	# these constants were used from SPIRE guide for 250, 350 and 500
	# and calculated using the equation that is explain in the end
	if wavelength == 70:
		data = data * 1./((np.pi*header['CDELT2']/180.)**2*1.e6)
		beam = 'nan' 
		sigma = all_data['stDev'].data 
		sigma = sigma * 1./((np.pi*header['CDELT2']/180.)**2*1.e6)
		coverage = all_data['coverage'].data

	elif wavelength == 160:
		data = data * 1./((np.pi*header['CDELT2']/180.)**2*1.e6) 
		data[data==0] = float('nan')
		beam = beamsize['160']
		sigma = all_data['stDev'].data 
		sigma = sigma * 1./ ((np.pi*header['CDELT2']/180.)**2*1.e6)
		coverage = all_data['coverage'].data*31.25 # 31.25 is used to convert the units to number of bolometers

	elif wavelength == 250:
		beam = beamsize['250']
		sigma = all_data['error'].data 
		coverage = all_data['coverage'].data 

	elif wavelength == 350:
		beam = beamsize['350']
		sigma = all_data['error'].data
		coverage = all_data['coverage'].data

	elif wavelength == 500:
		beam = beamsize['500']
		sigma = all_data['error'].data
		coverage = all_data['coverage'].data 	

	# change the wavelength and units in the header
	header['FIELD'] = (field, 'Name of the field')
	header['WAVELNTH'] = (wavelength, '[micrometer] The reference wavelength' )                           
	header['BUNIT'] = ('MJy/sr ','Unit of the data' )
	header['RESO'] = (beam, '[arcsec] Resolution of the map')
	header['POSANGLE'] = (header_true['POSANGLE'], '[deg] Position Angle of pointing')

	try:
		if header_true.comments['META_8'] == '[MJy sr-1] Offset added to SPIRE image as ob&':
			header['H_OFFSET'] = (header_true['META_8'], '[MJy/sr] offset, as derived by Herschel')
			header['H_OFF_E'] = (header_true['META_9'], '[MJy/sr] offset error, as derived by Herschel')
		elif header_true.comments['META_7'] == '[MJy sr-1] Offset added to SPIRE image as ob&':
			header['H_OFFSET'] = (header_true['META_7'], '[MJy/sr] offset, as derived by Herschel')
			header['H_OFF_E'] = (header_true['META_8'], '[MJy/sr] offset error, as derived by Herschel')
	except:
		print ('There is no offset provided by Herschel')

	try:
		header.remove(keyword='QTTY____')
		header.remove(keyword='INFO____')
	except:
		print()	

	# get the path and save the fits file
	if wavelength == 160:
		path  = file_name[:-5] 
	else:
		path  = file_name[:-8] 	

	path_image = '{0}/{1}-{2}.image.fits'.format(folder,field,int(wavelength))
	print (path_image)
	os.system('rm {0}'.format(path_image))
	header['EXTNAME'] = ('Image','name of this file')

	fits.writeto(path_image, data, header)

	path_sigma = '{0}/sigma_maps/{1}-{2}.sigma.fits'.format(folder,field,int(wavelength))
	print (path_sigma)
	os.system('rm {0}'.format(path_sigma))
	header['EXTNAME'] = ('Sigma','name of this file')

	fits.writeto(path_sigma, sigma, header)

	header['BUNIT'] = ('1 ','Unit of the data' )

	path_coverage = '{0}/coverage_maps/{1}-{2}.coverage.fits'.format(folder,field,int(wavelength))
	print (path_coverage)
	os.system('rm {0}'.format(path_coverage))
	header['EXTNAME'] = ('Coverage','name of this file')

	fits.writeto(path_coverage, coverage, header)
	return 


###################################################################################################
################################## MAIN PROGRAM STARTS HERE #######################################
###################################################################################################

main_folder = path_to_folder

# list of wavelengths
wavelength_list = [160. , 250., 350., 500.]
field = sys.argv[1] 	# Name of the field
print ('Field:', field)

folder = '{0}/herschel_data/{1}'.format(main_folder,field)

os.system('mkdir {0}/sigma_maps'.format(folder))
os.system('mkdir {0}/coverage_maps'.format(folder))
print ()
# changing units. 
all_files = sorted(glob.glob('{0}/*ext.fits'.format(folder)))
for j,i in enumerate(wavelength_list):
	print ('Wavelength:',i)
	file_name = all_files[j]
	unit_change(file_name, field, i, 'same')
	print ()

	
'''
to change from Jy/pixel to MJy/sr:

1 [sr] = ((180/pi)*3600)^2 [arcsec^2] = (C [arcsec^2])

1 [Jy/pixel] = 1 Jy/(pixelsize^2 [arcsec^2]) (pizelsize is 3" for 160 and 2" for 70)

1 [Jy/pixel] = 1 Jy/(pixelsize^2 [arcsec^2]) * (C [arcsec^2])/[sr] * [MJy]/10^6[Jy]

1 [Jy/pixel] = ((180/pi)*3600)^2 /(pixelsize^2 * 10^6) [My/sr]

'''

