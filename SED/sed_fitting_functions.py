from __future__ import print_function
###################################################################################################
'''
This script has functions that are used in program: 'creating_tau_and_temp.py' to generate tau at 1 THz,
column density of H2 and dust temperature maps. 

This script requires following functions and files:
initial_parameters.py
colour_correction.py
stored_offsets.py

The programs consists of following functions: 
- blackbody 			: normal blackbody
- modified_blackbody	: modified blackbody
- make_SED_maps			: Spectral Energy Distrubution fitting function (this uses peval and residuals)
- getChiSquare			: get Chi square value using covariance matrix 
- intensity      		: Calculate blackbody for least square funtion

For more information look in each function


By: Ayushi Singh
Last Modified: 24 April 2020
'''
###################################################################################################

# package required 
import numpy as np
from scipy.optimize import curve_fit
from astropy import constants as const
from astropy.io import fits

# other scripts 
from initial_parameter import field_name, nu0, kappa0, gastodust, muh2
from initial_parameter import temp_file, tau_file, beta_file, beta
from initial_parameter import other_uncertainty, types_of_errors, offset, astrometric_error
import colour_correction as cc
import stored_offsets as so

# Supress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)
###################################################################################################
#################################### CONSTANTS ####################################################
###################################################################################################


"""

http://arxiv.org/pdf/1101.4654v1.pdf

Constants
----------
muh2 : float
   	The mass per molecule of H2.  Usually is 2.8.       	  	
h : float
	Planck constant in erg.s
k : float
	Boltzmann constant in erg/K
c : float
	Speed of light in cm/s
mh : float
	Mass of hydrogen atom in g
mass_sun : float
	Mass of sun in g
pc : float 		 	
	conversion value of 1 pc to cm		
		
"""
           				

h,k,c,mh = const.h.cgs , const.k_B.cgs , const.c.cgs, const.m_p.cgs*1.00794
h,k,c,mh = h.value, k.value, c.value, mh.value	# [erg.s], [erg/K], [cm/s], [g]
mass_Sun, pc = const.M_sun.cgs, const.pc.cgs
mass_Sun, pc = mass_Sun.value, pc.value 		# [g], [cm]

Mjy = 1e-17	 	 								# divide this to convert [erg/sr.s.Hz.cm^2] to [MJy]

###################################################################################################
#################################### FUNCTIONS ####################################################
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

def modified_blackbody(nu, temperature, Tau_nu_0):
	"""
	Give the intensity value for a given frequency, temperature, log(N)
	using modified blackbody equation: Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1  (1 - e^(-tau_nu))
	To get the simple blackbody, it uses the blackbody function

	Parameters
	----------
	nu : float
		frequency in Hz
	temperature : float
		temperature in Kelvin
	Tau_nu_0 : float
		the optical depth		

	Returns
	-------
	I : float
	    Intensity by fitting modified blackbody in MJy/sr

	For other constants: look where they are defined at the top    
	"""
	#Tau_nu_0 =  muh2 * mh * (kappa0 / gastodust) *  10.0**logN

	tau =Tau_nu_0*(nu/nu0)**beta_val 

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

def run_SED_leastsq(wavelength, im, cpd, gradients, save_all=False):
	'''
	This function uses wavelengths and images at those wavelengths to 
	generate column density and temperature values at each pixel using
	Spatical Energy Distrubution fitting. This uses following functions:
	peval, residuals, modified_blackbody and blackbody.

	Parameters
	----------
	wavelengths : 1-D array
		array of all the wavelength in micron	
	im : 3-D array
		flux images in MJy/sr
		coorediantes are number of images corresponding to number of 
		wavelengths and two spatical coordinates 
	cpd : 3-D array
		instrument and gradient based uncertatinty maps 		
	gradients: 3-D array
		gradient maps of the data

	Returns
	-------
	Temperature : 2-D array
	    map of temerature in Kelvin
	tau_nu0: 2-D array
	    map of optical depth
	MassCD: 2-D array
	    map of column density in M_sun/pc^2    
	'''
	print ('RUNNING SED FITTING SCRIPT')

	im1 = im[0]

	if temp_file != None:
		planck_temp = fits.getdata(temp_file)
	
	if tau_file != None:		
		planck_tau = fits.getdata(tau_file)

	if beta_file != None:	
		planck_beta = fits.getdata(beta_file)
	
	# calculate frequency from each wavelength
	Wavelength = np.array(wavelength).astype('float64')  
	print ('Wavelength in [micron]:', Wavelength) 
	Wavelength = Wavelength/10000				# [cm]
	Frequency = (c/Wavelength)					# [Hz] 
	print ('Frequency in [GHz]:', Frequency/1e9)
	
	# temperature [K] : an empty map
	Temperature = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 

	# column-density [log[cm^-2]]: an empty map	
	tau_nu0 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 

	# chi square: an empty map	
	Chi_sq = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	
	
	# SED final: an empty map	
	sed160 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	sed250 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	sed350 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	sed500 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]])
	
	err160 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	err250 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	err350 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	err500 = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	

	# errors on temp and tau: empty maps
	T_err = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 
	tau0_err = np.zeros([np.shape(im1)[0],np.shape(im1)[1]]) 

	# this will count number of total pixel		
	pixel= 0						 
	
	#### for loop with i goes over the x values of the pixels and 
	#### for loop with j goes over the y values of the piexls and
	#### for loop with q goes over each of the 4 or n imgaes
	if offset == True:
		offset_value = so.getOffset(field_name)
	else:
		offset_value = np.zeros(len(Frequency))	

	print ('\nThe offset values subtracted from Data is: ', offset_value, '\n')
	for w in range(len(Frequency)):
		im[w]=im[w]-offset_value[w] 		# [MJy/sr]


	for i in range(np.shape(im1)[0]):
		print ("Running pixel row x: ", str(i+1), 'out of ', str(np.shape(im1)[0]), end="\r")
		for j in range(np.shape(im1)[1]):
			# will hold all the flux for that pixel
			Flux = np.zeros(len(Frequency)) 

			global beta_val
			if beta == None:
				beta_val = planck_beta[i][j]#gv.get_beta(i,j)
			else:
				beta_val = beta

			global Err
			Err = np.zeros(len(Frequency))	
			
			global Cov
			Cov = np.zeros([len(Frequency), len(Frequency)]) 

			'''
			if offset == True:
				offset_value = so.getoffset(field_name)
			else:
				offset_value = np.zeros(len(Frequency))	
			'''
			for q in range(len(Frequency)):
				Flux[q]=im[q][i][j]#-offset_value[q] 		# [MJy/sr]
			
			for q in range(len(Frequency)):
				if np.isnan(Flux[q]) == True or Flux[q]<=0.0:
					break
				elif np.isnan(cpd[q][i,j]) == True:
					break
				else:
					lmbda = str(int(wavelength[q]))

					if types_of_errors['coverageBased'] == True:
						Cov[q,q] = Cov[q,q]+cpd[q][i][j]**2

					if types_of_errors['undulation'] == True:
						Cov[q,q] = Cov[q,q]+(other_uncertainty[lmbda]['undulation'])**2
					    
					if types_of_errors['calibration'] == True:
						Cov[q,q] = Cov[q,q]+(other_uncertainty[lmbda]['calibration']*Flux[q])**2
					    
					if types_of_errors['corr_cali'] == True:
						for f in range(len(Frequency)):
							pacs_spire_per = 0.0
							spire_spire_per = 1.0
							factor = np.array([[ spire_spire_per,pacs_spire_per,pacs_spire_per,pacs_spire_per],
												[pacs_spire_per, spire_spire_per, spire_spire_per, spire_spire_per],
												[pacs_spire_per, spire_spire_per, spire_spire_per, spire_spire_per],
												[pacs_spire_per, spire_spire_per, spire_spire_per, spire_spire_per]])			
							lmbda2 = str(int(wavelength[f]))
							if q == f:
								Cov[q,f] = (Cov[q,f]+(other_uncertainty[lmbda]['corr_cali']*Flux[q])**2)
							else:
								sigma1 = (other_uncertainty[lmbda]['corr_cali']*Flux[q])
								sigma2 = (other_uncertainty[lmbda2]['corr_cali']*Flux[f])
								Cov[q,f] = Cov[q,f]+(factor[q,f]*sigma1*sigma2)  

					if types_of_errors['ciba'] == True:
						for f in range(len(Frequency)):
							
							factor = np.array([[1.00,0.95,0.86,0.80],
									   [0.95,1.00,0.95,0.86],
									   [0.86,0.95,1.00,0.95],
									   [0.80,0.86,0.95,1.00]])
							lmbda2 = str(int(wavelength[f]))
							if q == f:
								Cov[q,f] = Cov[q,f]+(other_uncertainty[lmbda]['ciba'])**2
							else:
								sigma1 = other_uncertainty[lmbda]['ciba']
								sigma2 = other_uncertainty[lmbda2]['ciba']
								Cov[q,f] = Cov[q,f]+(factor[q,f]*sigma1*sigma2)

					if types_of_errors['astrometric'] == True:
						for f in range(len(Frequency)):
							
							factor = np.array([[1.00,-1.00,-1.00,-1.00],
									   [-1.00,1.00,1.00,1.00],
									   [-1.00,1.00,1.00,1.00],
									   [-1.00,1.00,1.00,1.00]])
							lmbda2 = str(int(wavelength[f]))
							if q == f:
								Cov[q,f] = Cov[q,f]+(astrometric_error*gradients[q][i][j])**2
							else:
								sigma1 = astrometric_error*gradients[q][i][j]
								sigma2 = astrometric_error*gradients[f][i][j]
								Cov[q,f] = Cov[q,f]+(factor[q,f]*sigma1*sigma2)

				Err[q] = np.sqrt(Cov[q,q])

			Err[Err <=0.0] = float('nan')
			Flux[Flux <=0.0] = float('nan')	
			is_there_any_nan = np.isnan(Flux)
			
			if np.any(is_there_any_nan) == True or np.any(np.isnan(Err)) == True:
				# store values at corresponding pixel
				Temperature[i][j] = float('nan')	 

				tau_nu0[i][j] = float('nan')

				Chi_sq[i][j] = float('nan')

				if save_all == True:
				
					sed160[i][j] = float('nan')
					sed250[i][j] = float('nan')
					sed350[i][j] = float('nan')
					sed500[i][j] = float('nan')

					err160[i][j] = float('nan')
					err250[i][j] = float('nan')
					err350[i][j] = float('nan')
					err500[i][j] = float('nan')
				

				T_err[i][j] = float('nan')
				tau0_err[i][j] = float('nan')

			else:	
				x = Frequency 				# [Hz]
				y = Flux 					# [MJy/sr]

			
				# these are intial guess and can be changed at the top of the program with other constants. 
				# NOTE: in some cases these intitial guesses can make a big difference. 
				if temp_file == None:
					temp_guess = 17.
				else:
					temp_guess = planck_temp[i][j]

				if tau_file == None:	
					Tau_nu0_guess = 0.001
				else:
					Tau_nu0_guess = planck_tau[i][j]*(nu0/353.0e9)**beta_val	

				p_initial = [temp_guess , Tau_nu0_guess]	
				try:
					# run the curve_fit function that minimized the chi square: NOTE: peval is a function (see below) 
					p_final, p_error = curve_fit(intensity, x,y, p0=p_initial, sigma = Cov, bounds=(0, [100., 1.]), absolute_sigma = True)

					# calculate predicted intensity values
					y_final = intensity(x,*p_final)

					# getting errors on the parameters
					perrors = np.sqrt(np.array([p_error[0][0],p_error[1][1]]))

					# calcualte chi^2
					chi2 = getChiSquare(y, y_final, Cov)

					# store values at corresponding pixel
					Temperature[i][j] = p_final[0]	 

					tau_nu0[i][j] = p_final[1]

					Chi_sq[i][j] = chi2

					if save_all == True:

						sed500[i][j] = y_final[3]
						sed350[i][j] = y_final[2]
						sed250[i][j] = y_final[1]
						sed160[i][j] = y_final[0]			
						
						err500[i][j] = Err[3]
						err350[i][j] = Err[2]
						err250[i][j] = Err[1]
						err160[i][j] = Err[0]
					
					T_err[i][j] = perrors[0]
					tau0_err[i][j] = perrors[1]
					
				except:
					# store values at corresponding pixel
					Temperature[i][j] = float('nan')	 

					tau_nu0[i][j] = float('nan')

					Chi_sq[i][j] = float('nan')

					if save_all == True:
					
						sed160[i][j] = float('nan')
						sed250[i][j] = float('nan')
						sed350[i][j] = float('nan')
						sed500[i][j] = float('nan')

						err160[i][j] = float('nan')
						err250[i][j] = float('nan')
						err350[i][j] = float('nan')
						err500[i][j] = float('nan')
					

					T_err[i][j] = float('nan')
					tau0_err[i][j] = float('nan')

				'''
				print ()
				print ('~~~~~~~~~~~~~~~~~~~~~~~~~~')
				print ('Flux: ', Flux)
				print ('Err: ', Err)

				print ('Tau: ', tau_nu0[i][j])
				print ('Err Tau: ', tau0_err[i][j] )

				print ('Temp: ', Temperature[i][j])
				print ('Err Temp: ', T_err[i][j])

				print ('Chisq: ', Chi_sq[i][j])
				print ('~~~~~~~~~~~~~~~~~~~~~~~~~~')
				'''

				'''
				if pixel % 50000 ==0:# and plots == True:
					import matplotlib.pyplot as plt
					plt.ion()
					plt.figure()
					print ('pixel number:',pixel)
					freq = np.arange(2,4,0.01)
					freq=10**freq
					plt.loglog(x,y,'.', label='Data')
					plt.loglog(freq, intensity(freq,*p_final[0]), label='Blackbody Fit')
					plt.xlabel('Frequency [GHz]')
					plt.ylabel('Flux [Jy/sr]')
					plt.title('{0}/{1}'.format(p_final[0], pixel))
					plt.legend(loc=0)
					plt.show()
				'''
			pixel+=1	
	
	#### this creats a map of column density in the units of Msun/pc^2
	#tau_nu0 =  muh2 * mh * (kappa0 / gastodust) *  10.0**logN

	C_D = tau_nu0/(muh2 * mh * (kappa0 / gastodust)) # this is divided by (nu/nu0)^beta but because nu = nu0, this is just 1

	if save_all == True:
		return Temperature, tau_nu0, C_D, Chi_sq, T_err, tau0_err, sed160,sed250,sed350,sed500, err160, err250, err350, err500
	else:
		return Temperature, tau_nu0, C_D, Chi_sq, T_err, tau0_err

def intensity(x,temp, tau):
	'''
	function that calculate blackbody of list of intensities using
	modified_blackbody and blackbody functions after they have been colour corrected. 

	Parameters
	----------
	x : 1-D array
		frequency in Hz
	temp : float
		temperature in K
	tau : float
		optical depth

	Returns
	-------
	I : 1-D array
	    Intensity for each frequency by fitting modified blackbody in MJy/sr

	'''
	# get the colour correction values
	output = np.zeros(len(x))
	cc_vals = cc.colourcorrection(temp, beta_val,frequency=x)
	#cc_vals = output+1.

	# this loop goes through 
	for i in range(len(x)):
		predicted_val = (modified_blackbody(x[i], temp, tau))
		cc_predicted_val = predicted_val/cc_vals[i]
		output[i] = cc_predicted_val
	return output		


def getChiSquare(y, y_final, Cov):
	'''
	Chi Square function

	Parameters
	----------
	y : float
		intensity given in MJy/sr
	x_final : float
		model intensity given in MJy/sr
	Cov_inv: 1-D array		
		Inverse of Covariance Matrix.

	Returns
	-------
	Chisq : float
	    Chi square value
	'''
	try: 
		Cov_inv = np.linalg.inv(Cov)
		delta = y - y_final
		multi_a = np.matmul(delta,Cov_inv)		
		Chisq = np.matmul(multi_a, delta.T)

	except:
		Chisq = float('nan')

	#print(Chisq)
	return Chisq

