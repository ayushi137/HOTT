from __future__ import print_function
##############################################################################
'''
Python script to get model intensity map of cloud using Planck map of
Tau, Temperature and Beta for any given wavelength. 

This script requires following functions and files:
initial_parameter.py

To run:
run backward_sed_modelling.py 'Planck_Field'

Example:
run backward_sed_modelling.py TauLL

By: Ayushi Singh
Last Modified: 30 December 2019
'''
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import constants as const
import os, sys

# importing other scipts
from initial_parameter import path_to_folder
from initial_parameter import tau353GHz, temperature, beta

"""
N = 1e22 is the column density in cm^-2

Other constants
----------
nu0 : float
	The frequency at which the opacity power law is locked. In this case, it's 353 GHz. 
	This correspond to the frequency of tau0 planck maps.
h : float
	Planck constant in erg.s
k : float
	Boltzmann constant in erg/K
c : float
	Speed of light in cm/s
"""
                 
nu0= 353e9 								# [Hz] 
h,k,c = const.h.cgs , const.k_B.cgs , const.c.cgs 
h,k,c = h.value, k.value, c.value		# [erg.s], [erg/K], [cm/s]

Mjy = 1e-17	 	 						# divide this to convert [erg/sr.s.Hz.cm^2] to [MJy]

###################################################################################################
###################################################################################################
###################################################################################################

def blackbody(nu,temperature):
	"""
	Give the intensity value for a given frequency and temperature using 
	blackbody equation: Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1. 

	Parameters
	----------
	nu : float
		frequency in Hz
	temperature : float
		temperature in Kelvin

	Returns
	-------
	I : float
	    Intensity by fitting blackbody in erg/(sr.s.Hz.cm^2)

	"""
	# calculation is broken down in parts 
	aa = h*nu
	ab = (k*temperature)
	a = (np.exp(aa/ab) - 1)

	I = (2*h*nu**3 / c**2 )* a**-1

	return I 											# [erg/sr.s.Hz.cm^2]

############### function that calculate modified blackbody #############################

def modified_blackbody(nu, temperature, tau0, Beta):
	"""
	Give the intensity value for a given frequency, temperature, tau at 353 GHz and Beta
	using modified blackbody equation: Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1  (1 - e^(-tau_nu))
	To get the simple blackbody, it uses the function blackbody 

	Parameters
	----------
	nu : float
		frequency in Hz
	temperature : float
		temperature in Kelvin
	tau0 : float
		the optical depth at 353 GHz
	beta : float
		Dust emissivity: the spectral index that determines the opacity.		

	Returns
	-------
	I : float
	    Intensity by fitting modified blackbody in MJy/sr

	For other constants: look where they are defined at the top    
	"""
	tau = tau0 * (nu/nu0)**Beta 

	# this is the modification due to the optical depth 
	modification = (1.0 - np.exp(-1.0 * tau))

	# this calls blackbody function that just give plain blackbody intensity
	I = blackbody(nu, temperature)

	# intensity is modified 
	I = I*modification   								# [erg/sr.s.Hz.cm^2]

	# the I is divied by 1e-23 this changes [erg/(s.cm^2.Hz)] to Jy
	# the I is divied by 1e-17 this changes [erg/(s.cm^2.Hz)] to MJy
	I = I/Mjy										    # [MJy/sr]


	return I

def modelling (nu,tau0,T,Beta):
	"""
	Give the intensity value for a given frequency, temperature, tau at 353 GHz and Beta
	using modified blackbody equation: Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1  (1 - e^(-tau_nu))
	To get the simple blackbody, it uses the function blackbody 

	Parameters
	----------
	nu : float
		frequency in Hz
	temperature : float
		temperature in Kelvin
	tau0 : float
		the optical depth at 353 GHz
	beta : float
		Dust emissivity: the spectral index that determines the opacity.		

	Returns
	-------
	I : n-D array
	    Intensity by fitting modified blackbody in MJy/sr for each pixel of 2-D array 
	    and n numbers of wavelength    
	"""
	# total number of wavelengths
	num_wavelength = len(nu)

	# empty arrays for the intensity will store unit [MJy/sr] 
	# this has 3 dimension shape (number of wavelengts, length in x, length in y)
	Intensity = np.zeros([num_wavelength, np.shape(T)[0],np.shape(T)[1]])

	# for each pixel take hte temp, beta nad tau and apply blackbody fitting
	for i in range(np.shape(T)[0]):
		print ("Running pixel row x: ", str(i+1), 'out of ', str(np.shape(T)[0]), end="\r")
		for j in range(np.shape(T)[1]):
			temp = T[i][j]
			B = Beta[i][j]
			TAU0 = tau0[i][j]
			for k in range(len(nu)):
				Intensity[k][i][j] = modified_blackbody(nu[k], temp,TAU0,B)
	return Intensity

###################################################################################################
################################## MAIN PROGRAM STARTS HERE #######################################
###################################################################################################

main_folder = path_to_folder

# NEED TO MODIFY
name = sys.argv[1]#'TauL'
print (name)
folder = main_folder+'/Planck'
print (folder)
result_folder = main_folder+'/Models'#working_scripts/offset_v3/planck_models/'+name
print (result_folder)
#os.system('mkdir {0}'.format(result_folder))

#wavelength = np.arange(100, 801, 1)
wavelength = np.array([70,100,160,250,350,500])
wavelength_name =  []
for i in wavelength:
	if i < 100:
		wavelength_name.append('0'+str(i))
	elif i >= 100:
		wavelength_name.append(str(i))		


# This gets the file path and name
# NOTE: "name" defined earlier should be as it is given to the file 
tau0 = '{0}/{1}{2}'.format(folder,name, tau353GHz)
T = '{0}/{1}{2}'.format(folder,name, temperature)
Beta = '{0}/{1}{2}'.format(folder,name, beta)

# upload the tau, beta and temp files
tau0 = fits.getdata(tau0)								# tau0 at 353 GHz
T, head = fits.getdata(T,header = True)					# temperature [K]
Beta = fits.getdata(Beta)								# beta

# get the header of Temp and modify to make it a header for intensity maps
head['BUNIT'] = 'MJy/sr'
head['BTYPE'] = 'FLUX'


# convert wavelength [micron] to frequency [GHz]
wave = wavelength*1./10000								# [cm]
frequency = (c/wave)									# [Hz]

# Get the intensity using modelling function
Intensity = modelling(frequency,tau0,T,Beta)

print ("\n names of all the output files are as follows,")

########## SAVING ##################
# save them all with appropriate file name
for i in range(len(frequency)):
	print ("file #:", i+1 , "of", len(frequency))
	print ('\t{0}/{2}_HFI_model_{1}.fits'.format(result_folder,wavelength_name[i], name))
	os.system('rm {0}/{2}_HFI_model_{1}.fits'.format(result_folder,wavelength_name[i], name))
	fits.writeto('{0}/{2}_HFI_model_{1}.fits'.format(result_folder,wavelength_name[i], name),Intensity[i], head)


'''
########## PLOTTING ##################

# Choose test pixel 
x = 100
y = 100 

plt.ion()

# this plots the all 6 maps for 70, 160, 250, 350, 500
plt.figure(figsize=(15,6))
maxes = [100,700,1800, 1250,800,300]
inten = []
for i in range(len(frequency)):
	inten.append(Intensity[i][x][y])
	plt.subplot(2,3,i+1)
	plt.imshow(Intensity[i], origin = 'lower', vmin=0, vmax=maxes[i])
	plt.colorbar()
	plt.title('{0} micron'.format(wavelength_name[i]))

inten = np.array(inten)
freq = np.arange(2,4,0.01)
freq=10**(freq+9)

test_inten = np.zeros(len(freq))
for i in range(len(freq)):
	test_inten[i] = modified_blackbody(freq[i], T[x][y],tau0[x][y],Beta[x][y])

plt.figure()
plt.loglog(frequency, (inten), 'or')
plt.loglog(freq, (test_inten), 'b')
'''
