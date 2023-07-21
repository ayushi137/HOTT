from __future__ import print_function
##############################################################################
'''
This script is used to create the common science area mask for Herschel data 
using both PACS and SPIRE coverage maps. This make sure that the mask incorporate 
all four scans: normal and orthogonal scans for SPIRE and PACS.  

This script requires following functions and files:
initial_parameter.py
required_functions.py
boundary.py --- Note: this is where all the edge values are stored so it doesn't have 
be entered manually every single time. Once the value has been generated when the regions 
was first analyzed, they need to be stored in bounday.py for future auto-run. 

To run: 
run cropping_funtion.py 'Planck Field' 'Field' 

Example:
run cropping_funtion.py TauLL TauS2 

By: Ayushi Singh
Last Modified: 30 December 2019
'''
##############################################################################
# importing modules
from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys, os

# importing other scipts
from initial_parameter import path_to_folder
import required_functions as rf
import boundary

# alternative to tight layout
rcParams.update({'figure.autolayout': True})
plt.ion()

# Supress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=UserWarning) 

# starting name of planck file
planck_name = sys.argv[1]

# starting name of herschel file
name = sys.argv[2]

print (name + '\n')

#####################################################################################################
######################################## FUNCTIONS ##################################################
#####################################################################################################
def linefit (a, b, x):
	m = np.float(b[1] - a[1])/np.float(b[0] - a[0])
	b = b[1] - (b[0]*m)
	y = x*m + b
	return int(y)

def cropping(cropped, top, bottom, left, right):

	for i in range(np.shape(cropped)[1]):
		t = linefit(top, left, i)
		l = linefit(left, bottom, i)
		b = linefit(bottom, right, i)
		r = linefit(right, top, i)
		for j in range(np.shape(cropped)[0]):
			if j >= t:
				cropped[j][i] = float('nan')
			if j <= l:	
				cropped[j][i] = float('nan')
			if j <= b:
				cropped[j][i] = float('nan')
			if j >= r:
				cropped[j][i] = float('nan')
			if np.isnan(cropped[j][i]) == False:
				cropped[j][i] = 1
	return cropped

#####################################################################################################
################################## VARIOUS FILE NAMES ###############################################
#####################################################################################################

main_folder = path_to_folder

# input file: coverage maps and any planck map with appropriate pixelsize 
filepacs = '{0}/herschel_data/{1}/coverage_maps/{1}-160.coverage.fits'.format(main_folder, name)
filespire_250 = '{0}/herschel_data/{1}/coverage_maps/{1}-250.coverage.fits'.format(main_folder, name)
filespire_350 = '{0}/herschel_data/{1}/coverage_maps/{1}-350.coverage.fits'.format(main_folder, name)	
filespire_500 = '{0}/herschel_data/{1}/coverage_maps/{1}-500.coverage.fits'.format(main_folder, name)
planck_file = '{0}/Models/{1}_HFI_model_160.fits'.format(main_folder,planck_name)

# ouput of coverage maps and regridded files
overlap = '{0}/coverage/{1}_overlap.png'.format(main_folder, name)

# output mask files
grid_planck = '{0}/grid/{1}_planck.fits'.format(main_folder, name)
grid_500 = '{0}/grid/{1}_500.fits'.format(main_folder, name)
grid_350 = '{0}/grid/{1}_350.fits'.format(main_folder, name)
grid_250 = '{0}/grid/{1}_250.fits'.format(main_folder, name)
grid_160 = '{0}/grid/{1}_160.fits'.format(main_folder, name)
grid = '{0}/grid/{1}_grid.png'.format(main_folder, name)


#####################################################################################################
################################## MAIN PROGRAM #####################################################
#####################################################################################################

# regriding spire_500 map to pacs
os.system("rm {0}".format(filespire_500[:-5]+'.resamp.fits'))
rf.regrid(filespire_500,filepacs, resultimage = 'same', header = None)

# get the new coverage map
pacs = fits.getdata(filepacs)
spire = fits.getdata(filespire_500[:-5]+'.resamp.fits')

# The SPIRE coverage maps provides with number of bolobeter per pixel. However, PACS is 
# in the units of pixel that is the size of the bolometer in the focal plane to get them 
# in the same unit, we need to multiply by 31.5

# modifying to make maps readable
spire[np.isnan(spire)] = 0
pacs[np.isnan(pacs)] = 0     
maxs = np.max(spire)
maxp = np.max(pacs)
pacs[pacs <= 0.0001] = float('nan')
spire[spire <= 0.0001] = float('nan')

# combining both images
combined = (spire/maxs) - (pacs/maxp)


# ploting the image that is used to create the mask
plt.figure()
plt.imshow(combined, origin='lower',vmin = -0.6, vmax = 0.6)
plt.grid(which = 'major')
plt.colorbar()

# saving the previous image as png
plt.savefig(overlap)

print ()

try:
	top,left,bottom,right = boundary.boundaries(name, 'edge')
except UnboundLocalError:
	print ('The edge value are not stored in boundary.py. Enter them manually.')
	print()

	print ("If you see the image. Click on the image to see it's interactive.")

	now = input('Then press enter.')

	print ()
	see = input('It the image interactive? y or n: ')

	if see == 'n':
		print ("Run the field again individually in ipython. Don't run them using run_offset.py")
		print()
		sys.exit("The code has been terminated.")
	elif see == 'y':	
		print ("It's a diamond shape. Start with the top point and then go counter clockwise.")
		tx = input('top-right point x: ')
		ty = input('top-right point y: ')
		lx = input('left-top point x: ')
		ly = input('left-top point y: ')
		bx = input('bottom-left point x: ')
		by = input('bottom-left point y: ')
		rx = input('right-bottom point x: ')
		ry = input('right-bottom point y: ')


		top = [tx,ty]
		left = [lx, ly]
		bottom = [bx,by]
		right = [rx,ry]

		print ()
		print ('Make sure to save it in the bounday.py for future reference. With the name "'+name+'" and the Purpose "Edge".')
	else:
		print ("Incorrect imput. Only say 'y' for yes or 'n' for no.")
		print()
		sys.exit("The code has been terminated.")



print ()
print ('\t'+'top    =',top ) 
print ('\t'+'left   =',left )
print ('\t'+'bottom =',bottom )
print ('\t'+'right  =',right )

# this calls the funtion that does the cropping using the values given above
cropped = np.zeros(np.shape(combined))
cropped = cropping(cropped, top, bottom, left, right)

# plot the cropped image
plt.figure()
plt.imshow(cropped, origin='lower')
plt.colorbar()

plt.savefig(grid)

#####################################################################################################
################################## SAVING OUTPUT ####################################################
#####################################################################################################

headerp = fits.getheader(filepacs)
headerp['EXTNAME'] = ('CSR Grid','Grid for Common Science Region')
headerp['BUNIT'] = ('None','Unit of the data' )

try:
	headerp.remove(keyword='WAVELNTH')
	headerp.remove(keyword='RESO')
	headerp.remove(keyword='H_OFFSET')
	headerp.remove(keyword='H_OFF_E')
except:
	print ()

# save the croppped image for 160 resoltion
os.system("rm {0}".format(grid_160))
fits.writeto(grid_160, cropped, headerp)

# regriding the cropped image to 500, 350, 250 and planck resolution
os.system("rm {0}".format(grid_500))
rf.regrid(grid_160,filespire_500, resultimage = grid_500, header = None)

os.system("rm {0}".format(grid_350))
rf.regrid(grid_160,filespire_350, resultimage = grid_350, header = None)

os.system("rm {0}".format(grid_250))
rf.regrid(grid_160,filespire_250, resultimage = grid_250, header = None)

os.system("rm {0}".format(grid_planck))
rf.regrid(grid_160,planck_file, resultimage = grid_planck, header = None)

