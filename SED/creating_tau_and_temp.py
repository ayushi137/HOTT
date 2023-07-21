from __future__ import print_function
###################################################################################################
'''
This program take the Herschel map that have been corrected for offsets and calculates 
the optical depth, column denstiy and dust temperature maps. The data set are of 160, 250, 350 and 
500 micron flux maps and their respesctive coverage maps.

This script requires following functions and files:
initial_parameters.py
prepare_maps.py
sed_fitting_functions.py

To run: creating_tau_and_temp.py Field
 
Example: creating_tau_and_temp.py TauS2

By: Ayushi Singh
Last Modified: 31 December 2019
'''
###################################################################################################

# package required 
import os, sys, time, datetime
import numpy as np
from astropy.io import fits
from astropy import constants as const

# other scripts 
from initial_parameter import field_name, folder, result_folder
from initial_parameter import wavelength, data_files, coverage_files
#from initial_parameter import file160, file250, file350, file500
#from initial_parameter import coverage160, coverage250, coverage350, coverage500
from initial_parameter import convolve_size, regrid_size, beamsize, types_of_errors, offset
from initial_parameter import nu0, kappa0, gastodust, muh2, beta, beta_file

import prepare_maps as pm
import sed_fitting_functions as sff

# Supress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

###################################################################################################

start = time.time()

print ('Running the SED fitting pipeline')

print ('\nField:', field_name)
print ()

os.system('mkdir '+result_folder)
print ('\nMade resulting Folder:', result_folder)

#wavelength = [160, 250, 350, 500]
#data_files = [file160, file250, file350, file500]	
#coverage_files = [coverage160, coverage250, coverage350, coverage500]	

print ('\nWavelength:', wavelength)

print ('\nFiles:', data_files)

header_file = data_files[np.argwhere (np.array(wavelength) == regrid_size )[0][0]]
print ('\nHeader File:', header_file)

# Convolving and regridding intensity maps
data_file_name = pm.generateMaps(wavelength, data_files, running=True)
print ()

# Convolving and regridding coverage maps
coverage_file_name = pm.generateMaps(wavelength, coverage_files, running=True)
print ()

# Calculating Gradient maps
grad_file_name = pm.generateGradientMaps(data_file_name, running=True)
print ()

# Calcualting Coverage Predicted Dispersion
uncertainty_file_name = pm.generateCPDMaps(wavelength ,data_file_name, coverage_file_name, grad_file_name, running=True)
print ()

# loading all the intensity and CPD maps
images = []
cpd = []
gradients = []
for i in range(len(data_file_name)):
	#print (data_file_name)
	#print (uncertainty_file_name)
	images.append(fits.getdata(data_file_name[i]))
	cpd.append(fits.getdata(uncertainty_file_name[i]))
	gradients.append(fits.getdata(grad_file_name[i]))

# Running the SED fitter
print ('\nWill save the results in:', result_folder)
print ('Applying the following errors:')
print ('- Coverage Predicted Dispersion:\t', types_of_errors['coverageBased'])
print ('- Undulations:\t\t\t\t', types_of_errors['undulation'])
print ('- Calibrations:\t\t\t\t', types_of_errors['calibration'])
print ('- Correlated Calibrations:\t\t', types_of_errors['corr_cali'])
print ('- CIBA:\t\t\t\t\t', types_of_errors['ciba'])
print ('- Astrometric:\t\t\t\t', types_of_errors['astrometric'])
print ('- Offset added:\t\t\t\t', offset)
print ()


### running data
get_all_the_maps = True
all_results = sff.run_SED_leastsq(wavelength, images, cpd, gradients, save_all= get_all_the_maps)
if get_all_the_maps == True:
	Temperature, tau_nu0, C_D, Chi_sq, T_err, tau0_err, sed160,sed250,sed350,sed500, err160, err250, err350, err500 = all_results
else:
	Temperature, tau_nu0, C_D, Chi_sq, T_err, tau0_err= all_results

# fixing up the header 
fitheader = fits.getheader(header_file)
#fitheader = fits.getheader(locals()["file"+str(regrid_size)])
try:
	fitheader.remove(keyword='WAVELNTH')
	fitheader.remove(keyword='OFFSET')
	fitheader.remove(keyword='OFFSET_E')
	fitheader.remove(keyword='SCALE')
	fitheader.remove(keyword='SCALE_E')

except:	
	print ('Zero-point Header remove failed')

try:
	fitheader.remove(keyword='PLANEX') 
	fitheader.remove(keyword='PLANEX_E') 
	fitheader.remove(keyword='PLANEY') 
	fitheader.remove(keyword='PLANEY_E') 
	fitheader.remove(keyword='PLN_CR1') 
	fitheader.remove(keyword='PLN_CR2') 

except:	
	print ('Plane Header remove failed')

try:
	fitheader.remove(keyword='H_OFFSET')
	fitheader.remove(keyword='H_OFF_E') 
	
except:	
	print ('SPIRE offset Header remove failed')

try: 
	fitheader.remove(keyword='DATA____')
	fitheader.remove(keyword='CLASS___')
except:	
	print ()

try:
	for i in range(7):
		fitheader.remove(keyword='        ') 
except:	
	print ()


fitheader['RESO'] = beamsize[str(regrid_size)]
if beta == None:
	fitheader['BETA'] = ('from Planck Model', 'Power law index of dust emissivity')
else:
	fitheader['BETA'] = (beta, 'Power law index of dust emissivity')
fitheader['AUTHOR'] = ('Ayushi Singh and Peter Martin')
fitheader['DATE'] = (str(datetime.date.today()), 'Date Created')
fitheader['ARTICLE'] = ('Singh & Martin (2022)', 'Reference publication')
fitheader['PRODUCT'] = ('HOTT maps', 'www.cita.utoronto.ca/HOTT')

print ('\n\nSaving Files')

# keeping only column density that is between 0 and 10^24 cm^-2
C_D[np.isnan(C_D)] = 0
C_D[np.isinf(C_D)] = 0
C_D[C_D<0]=float('Nan')
C_D[C_D>5e26]=float('Nan')

tau_nu0[np.isnan(C_D)] = 0
tau_nu0[np.isinf(C_D)] = 0

Temperature[np.isnan(Temperature)] = 0
Temperature[np.isinf(Temperature)] = 0
Temperature[Temperature<0]=float('Nan')
Temperature[Temperature>=200]=float('Nan')

beta_map = fits.getdata(beta_file)
Radiance = pm.create_radiance(Temperature, tau_nu0,beta_map ,nu0)

fitheader['BUNIT'] = 'K'
fitheader['HEADER'] = 'Dust Temperature'
os.system('rm '+'{0}/{1}_hott_temperature_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_temperature_orig.fits'.format(result_folder, field_name),Temperature,fitheader)

fitheader['HEADER'] = 'TAU at 1 THz'
os.system('rm '+'{0}/{1}_hott_tau1thz_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_tau1thz_orig.fits'.format(result_folder, field_name),tau_nu0,fitheader)

fitheader['BUNIT'] = 'None'
fitheader['HEADER'] = 'Chi_sqr'
os.system('rm '+'{0}/{1}_hott_chisquare_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_chisquare_orig.fits'.format(result_folder, field_name),Chi_sq,fitheader)

fitheader['BUNIT'] = 'K'
fitheader['HEADER'] = 'Dust Temperature error'
os.system('rm '+'{0}/{1}_hott_temperature_error_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_temperature_error_orig.fits'.format(result_folder, field_name),T_err,fitheader)

fitheader['BUNIT'] = 'None'
fitheader['HEADER'] = 'TAU at 1 THz error'
os.system('rm '+'{0}/{1}_hott_tau1thz_error_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_tau1thz_error_orig.fits'.format(result_folder, field_name),tau0_err,fitheader)

fitheader['BUNIT'] = 'Wm^-2sr^-1'
fitheader['HEADER'] = 'Radiance'
os.system('rm '+'{0}/{1}_hott_radiance_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_radiance_orig.fits'.format(result_folder, field_name),Radiance,fitheader)


fitheader['BUNIT'] = 'cm^-2'
fitheader['HEADER'] = 'H_2 Column Density'
fitheader['FREQ_0'] = (nu0/1e9, '[GHz] opacity power law lock frequency')
fitheader['KAPPA_0'] = (kappa0, '[cm^2/g] dust opacity at FREQ_0')
fitheader['GAS2DUST'] = (gastodust, 'gas to dust ratio')
fitheader['MU_H2'] = (muh2, '[amu] mean molecular mass of H2')
os.system('rm '+'{0}/{1}_hott_columndensity_orig.fits'.format(result_folder, field_name))
fits.writeto('{0}/{1}_hott_columndensity_orig.fits'.format(result_folder, field_name),C_D,fitheader)

os.system('rm '+'{0}/{1}_planck_beta_orig.fits'.format(result_folder, field_name))
os.system('cp {2} {0}/{1}_planck_beta_orig.fits'.format(result_folder, field_name, beta_file ))

fitheader.remove(keyword='FREQ_0')
fitheader.remove(keyword='KAPPA_0')
fitheader.remove(keyword='GAS2DUST')
fitheader.remove(keyword='MU_H2')
fitheader['BUNIT'] = 'MJy sr^-2'

if get_all_the_maps == True:
	fitheader['HEADER'] = 'model 160'
	os.system('rm '+'{0}/{1}_hott_model_160_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_model_160_orig.fits'.format(result_folder, field_name),sed160,fitheader)

	fitheader['HEADER'] = 'model 250'
	os.system('rm '+'{0}/{1}_hott_model_250_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_model_250_orig.fits'.format(result_folder, field_name),sed250,fitheader)

	fitheader['HEADER'] = 'model 350'
	os.system('rm '+'{0}/{1}_hott_model_350_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_model_350_orig.fits'.format(result_folder, field_name),sed350,fitheader)

	fitheader['HEADER'] = 'model 500'
	os.system('rm '+'{0}/{1}_hott_model_500_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_model_500_orig.fits'.format(result_folder, field_name),sed500,fitheader)

	fitheader['HEADER'] = 'error 160'
	os.system('rm '+'{0}/{1}_hott_error_160_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_error_160_orig.fits'.format(result_folder, field_name),err160,fitheader)

	fitheader['HEADER'] = 'error 250'
	os.system('rm '+'{0}/{1}_hott_error_250_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_error_250_orig.fits'.format(result_folder, field_name),err250,fitheader)

	fitheader['HEADER'] = 'error 350'
	os.system('rm '+'{0}/{1}_hott_error_350_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_error_350_orig.fits'.format(result_folder, field_name),err350,fitheader)

	fitheader['HEADER'] = 'error 500'
	os.system('rm '+'{0}/{1}_hott_error_500_orig.fits'.format(result_folder, field_name))
	fits.writeto('{0}/{1}_hott_error_500_orig.fits'.format(result_folder, field_name),err500,fitheader)


##### MASKING ALL THE MAPS TO CSR
print ('\n\nMasking for Common Science Region')
pm.generateMaskedMaps(result_folder)


print ("\n************ This field name was:", field_name, "*********************")
totalruntime = str(datetime.timedelta(seconds=time.time()-start))
print ('\nRun Time: ', str(datetime.timedelta(seconds=time.time()-start)))


try:
	os.system('echo SED finished for the field {0} in {1}. The results are saved in: {2} | mailx -r asingh@cita.utoronto.ca -s "Code Finished" ayushi.singh@mail.utoronto.ca'.format(field_name, totalruntime, result_folder))
except:
	pass
