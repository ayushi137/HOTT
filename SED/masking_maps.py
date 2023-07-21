from __future__ import print_function
###############
'''
Mask the data and save it in the masked subfolder

To run:
run masking_maps.py 'Field'

Example:
run masking_maps.py TauS2

By: Ayushi Singh
Last Modified: 7 August 2015
'''
###############

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys

plt.ion()

field = sys.argv[1]#'cepL1241'
print ("field:", field,"\n")
wavelength = '500'
path_to_mask = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project/grid'
mask = '{0}/{1}_{2}.fits'.format(path_to_mask, field, wavelength)

path_to_map = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project/tau_temp/{0}_wt_plane_10_resi_0'.format(field)
maps=glob.glob('{0}/*{1}*.fits'.format(path_to_map, wavelength))

os.system('mkdir {0}/{1}_masked'.format(path_to_map,field))

mask_data = fits.getdata(mask)
mask_data[mask_data<=0.001] = float('nan')
mask_data[mask_data>0.001] = 1

for i in maps:
	data, head = fits.getdata(i, header=True)
	data_masked = data*mask_data
	fits.writeto(i[:-10]+'.fits', data_masked, head)
	os.system('mv {0} {1}/{2}_masked'.format(i[:-10]+'.fits',path_to_map,field))

