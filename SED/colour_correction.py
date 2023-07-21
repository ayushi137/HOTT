##############################################################################
'''
This script has a function used to creating the colour correction interpolarator 
that will give the colour correction value for the given temperature and beta. 

This script requires following folder:
colour_correction: this folder has all the required files

By: Ayushi Singh
Last Modified: 4 January 2022 
'''
##############################################################################

# modules needed 
from astropy.io import fits
import numpy as np 
from scipy.interpolate import interp2d
from astropy import constants as const

# speed of light 
c = const.c.cgs
c = c.value # [cm/s]

# load all the required tables 
temp_160 = np.loadtxt('colour_correction/temperature_160.txt')
beta_160 = np.loadtxt('colour_correction/beta_160.txt')
temp = np.loadtxt('colour_correction/temperature.txt')
beta = np.loadtxt('colour_correction/beta.txt')
cc_160 = np.loadtxt('colour_correction/160_cc.txt')
cc_250 = np.loadtxt('colour_correction/250_cc.txt')
cc_350 = np.loadtxt('colour_correction/350_cc.txt')
cc_500 = np.loadtxt('colour_correction/500_cc.txt')

# interpolate each of the four tables
inter_160 = interp2d(temp_160,beta_160, cc_160, kind = 'cubic')
inter_250 = interp2d(temp, beta, cc_250, kind = 'cubic')
inter_350 = interp2d(temp, beta, cc_350, kind = 'cubic')
inter_500 = interp2d(temp, beta, cc_500, kind = 'cubic')

#####################################################################################################
######################################## FUNCTIONS ##################################################
#####################################################################################################

def colourcorrection(T, B,	wavelength=None, frequency=None):
	"""
	This gives a colour corrected value for a given temperature and beta.

	Parameters
	----------
	T : float
		Temperature in Kelvin
	B : float
		Beta, unit less
	wavelength : (optional) float or array or list
		wavelength in cm. If None then use frequency. Default is None.
	frequency : (optional) float or array or list
		frequency in Hz. If None then use wavelength. Default is None.	

	If both frequency and wavelength are None then it will 
	calulate for all 4 wavelength and returns the array. 		

	Returns
	-------
	cc : float or array or list
	    colour correction for each each of the four wavelength that the predicted value 
	    needs to be divided by to get the value what herschel would see.

	"""

	# if both wavelength and frequency is not given then give 
	# values for all 4 wavelengths
	if np.any(wavelength)==None and np.any(frequency)==None:

		cc = np.zeros(4)
		cc[0] = inter_160(T,B)
		cc[1] = inter_250(T,B)
		cc[2] = inter_350(T,B)
		cc[3] = inter_500(T,B)

	else:	
		# if frequency is given then first convert to wavelength
		if np.any(frequency) != None:
			frequency = np.array(frequency)
			wavelength = (c/frequency)					# [m]
			wavelength = wavelength*10000				# [cm]
		wave = wavelength

		# if wavelength an int or a float then run this part
		if  isinstance(wave,(int,float)):
			if wave >= 158. and wave <= 162.:
				cc = 1./inter_160(T,B)	
			elif wave >= 248. and wave <= 252.:	
				cc = inter_250(T,B)
			elif wave >=  348. and wave <= 352.:	
				cc = inter_350(T,B)
			elif wave >= 498. and wave <= 502.:	
				cc = inter_500(T,B)

		# if wavelength an list or array then run this part		
		else:
			length = len(wave)
			cc = np.zeros(length)
			for i in range(length):
				if wave[i] >= 158. and wave[i] <= 162.:
					cc[i] = 1./inter_160(T,B)	
				elif wave[i] >= 248. and wave[i] <= 252.:	
					cc[i] = inter_250(T,B)
				elif wave[i] >=  348. and wave[i] <= 352.:	
					cc[i] = inter_350(T,B)
				elif wave[i] >= 498. and wave[i] <= 502.:	
					cc[i] = inter_500(T,B)			

	return cc

'''
# test script
print colourcorrection (10,1.7)
print colourcorrection (10,1.7,wavelength=160.)	
print colourcorrection (10,1.7,frequency =1873702862500.0)

temp = np.arange(10,110,10)
for i in temp:
	print colourcorrection (i,1.71,wavelength=[160.,250.,350.,500.])	
#print colourcorrection (10,1.7,frequency =[  1.87370286e+12,   1.19916983e+12,   8.56549880e+11])
'''
	
