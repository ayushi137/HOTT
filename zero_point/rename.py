from __future__ import print_function
##############################################################################
'''
This script takes the raw file from the 'herschel_archive' folder and 
rename them and move them to appropriate folder called 'herschel_data' under
proper name. 
NOTE: Read the script before running it. Might need to change folder_pacs 
and/or folder_spire. Also make sure you unzip the .fits.gz files. 

To run:
run rename.py 'field' 'pacs_folder' 'spire_folder' ObservationID
Example:
run rename.py TauS2 Taurus/S2_p Taurus/S2_s 13 --- NOTE: the pac_folder is path to 
pacs data in herschel_archive folder. It's the folder that you get from Herschel 
archive.

By: Ayushi Singh
Last Modified: 30 December 2019
'''
##############################################################################
import numpy as np
import glob
import os, sys, gzip, shutil    

from initial_parameter import path_to_folder

main_folder = path_to_folder

field = sys.argv[1]   		# name of the field
pacs_folder = sys.argv[2]	# name of PACS folder in 'herschel_archive'
spire_folder = sys.argv[3]	# name of Spire folder in 'herschel_archive'
number = sys.argv[4]		# Observation ID

archive = 'herschel_archive'
data = 'herschel_data'

folder_pacs = '{0}/{1}/{2}/{3}*/level2_5/HPPJSM*R/'.format(main_folder,archive, pacs_folder,int(number))
folder_spire = '{0}/{1}/{2}/{3}*/level2_5/extd*W/'.format(main_folder,archive, spire_folder,int(number))

print (folder_pacs)
print (sorted(glob.glob(folder_pacs+'/*.fits.gz')))

result_folder = '{0}/{1}/{2}'.format(main_folder,data,field)

fits_pacs_gz = sorted(glob.glob(folder_pacs+'/*.fits.gz'))
fits_spire_gz = sorted(glob.glob(folder_spire+'/*.fits.gz'))

print ('All PACS .gz files:',fits_pacs_gz)
print ()
print ('All SPIRE .gz files:', fits_spire_gz)
print ()

# extracting all the files
try:
	with gzip.open(fits_pacs_gz[0], 'rb') as f_in:
		with open(fits_pacs_gz[0][:-3], 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
except:    
	print (fits_pacs_gz[0], ' does not exit.')

try:
	with gzip.open(fits_spire_gz[0], 'rb') as f_in:
		with open(fits_spire_gz[0][:-3], 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
except:    
	print (fits_spire_gz[0], ' does not exit.')

try:
	with gzip.open(fits_spire_gz[1], 'rb') as f_in:
		with open(fits_spire_gz[1][:-3], 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
except:    
	print (fits_spire_gz[1], ' does not exit.')

try:
	with gzip.open(fits_spire_gz[2], 'rb') as f_in:
		with open(fits_spire_gz[2][:-3], 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
except:    
	print (fits_spire_gz[2], ' does not exit.')

# getting the fits file
fits_pacs = sorted(glob.glob('{0}/*.fits'.format(folder_pacs)))
fits_spire = sorted(glob.glob('{0}/*.fits'.format(folder_spire)))

print ('All PACS:', fits_pacs)
print ()
print ('All SPIRE:', fits_spire)
print ()

# creating data folder 
print ('data folder:', result_folder)
os.system('mkdir {0}'.format(result_folder))


# renaming all the files

os.system('mv {0} {1}/{2}-160.notext.fits'.format(fits_pacs[0],result_folder, field))
os.system('mv {0} {1}/{2}-250.ext.fits'.format(fits_spire[2],result_folder, field))
os.system('mv {0} {1}/{2}-350.ext.fits'.format(fits_spire[1],result_folder, field))
os.system('mv {0} {1}/{2}-500.ext.fits'.format(fits_spire[0],result_folder, field))

