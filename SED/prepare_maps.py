from __future__ import print_function
###################################################################################################
'''
This script has functions that are used in program: 'creating_tau_and_temp.py' to generate more maps 
needed for the SED fitting. This includes the convolve and regdrided version of frequency and coverage 
maps, uncertainty maps and gradient maps. It also create the masked maps needed at the very end. 

This script requires following functions and files:
initial_parameters.py
required_functions.py

The programs consists of following functions:
- generateMaps 			: Get the convovled and regridded version of maps 
- generateGradientMaps	: Get the gradient maps 
- generateCPDMaps 		: Get the coverage predicted dispersion maps 
- generateMaskedMaps 	: Apply common science region to all the final maps

By: Ayushi Singh
Last Modified: 18 March 2022
'''
###################################################################################################
# package required 
import os, glob
import numpy as np
from astropy.io import fits
from astropy import constants as const
import scipy.signal as sig
import scipy.special as sp


# other scripts 
from initial_parameter import field_name, folder, path_to_folder, beamsize
from initial_parameter import convolve_size, regrid_size, to_get_errors, RIF_160, WN_160
import required_functions as rf

# Supress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 


def generateMaps(wavelength, files, running = True):
	'''
	This function uses wavelengths and images at those wavelengths to
	get convolved and regrided images of respective map.
	
	Parameters
	----------
	wavelengths : 1-D array
		array of all the wavelength in micron	
	file: list of strings
		Path and name of all the required files
	running: boolean
		If False, it will not create the new files but only 
		provide the name of the new files 

	Returns
	-------
	Final_names : List
		name of all the final files 

	'''
	print ()
	print ("CONVOLVING AND REGRIDDING \nTotal number of files:", len(files))

	# resolution for each wavelength
	original_res = []
	for i in wavelength:
		original_res.append(beamsize[str(i)])

	
	# get the template wavelength and file name
	# template is the one that is not changed in the end
	convolve_wavelength = convolve_size
	regrid_wavelength =  regrid_size
	template_map_name = files[np.argwhere(np.array(wavelength)==regrid_wavelength)[0][0]]
	required_res = beamsize[str(convolve_wavelength)]
	print ("Convolve template wavelength:", convolve_wavelength)
	print ("Convolve template resolution:", required_res)
	print ("Regrid template wavelength:", regrid_wavelength)
	print ("Regrid template image name:", template_map_name)
	print ()

	Final_names = []
	# this loops goes though each file 
	for i in range(len(files)):
		
		print ('Wavelength:',wavelength[i],'micron')

		# if the file is not 250 file then 
		# 1) it is convolved to the resolution of template map
		# 2) it is regrid to the shape and size of template
		# 3) saved into an array
		# else if, the file is the 250 file then 
		# saved into an array without alteration
		if wavelength[i] != convolve_wavelength and wavelength[i] != regrid_wavelength: 
			print ("For this file, going though convolving and regriding.")
			# file names used 
			conv_file = files[i][:-5]+'.conv.fits'			# convolve file
			resamp_file = conv_file[:-5]+'.resamp.fits'		# resampled image

			print (conv_file)
			print (resamp_file)
			
			if running == True:
				# remove if these files are already there
				os.system("rm {0}".format(conv_file))
				os.system("rm {0}".format(resamp_file))
				# convolve the image 
				rf.resconvolve(files[i],original_res[i],required_res, resultimage = conv_file)
			
				# make header
				head_orig = fits.getheader(conv_file)
				head_template = fits.getheader(template_map_name)

				try:
					head_template['OFFSET'] = head_orig['OFFSET']
					head_template['OFFSET_E'] = head_orig['OFFSET_E']
					head_template['SCALE'] = head_orig['SCALE']
					head_template['SCALE_E'] = head_orig['SCALE_E']
				except:
					pass

				try:
					head_template['H_OFFSET'] = head_orig['H_OFFSET']
					head_template['H_OFF_E'] = head_orig['H_OFF_E']
				except:
					pass
				
				# resample image 
				rf.regrid(conv_file, template_map_name, resultimage = resamp_file, header = head_template)
			
			Final_names.append(resamp_file)
			print ()	
		
		elif wavelength[i] == convolve_wavelength and wavelength[i] == regrid_wavelength: 
			print ("For this file, doing nothing.")

			Final_names.append(files[i])
			print ()

		elif  wavelength[i] == regrid_wavelength: 
			print ("For this file, going though only convolving.")
			# file names used 
			conv_file = files[i][:-5]+'.conv.fits'			# convolve file

			print (conv_file)
			
			if running == True:
				# remove if these files are already there
				os.system("rm {0}".format(conv_file))
				# convolve the image 
				rf.resconvolve(files[i],original_res[i],required_res, resultimage = conv_file)
			
			Final_names.append(conv_file)
			print ()	
	
		elif wavelength[i] == convolve_wavelength:
			print ("For this file, going though only regriding.")

			resamp_file = files[i][:-5]+'.resamp.fits'		# resampled image
			
			print (resamp_file)
			
			if running == True:
				os.system("rm {0}".format(resamp_file))
				# make header
				head_orig = fits.getheader(files[i])
				head_template = fits.getheader(template_map_name)

				try:
					head_template['OFFSET'] = head_orig['OFFSET']
					head_template['OFFSET_E'] = head_orig['OFFSET_E']
					head_template['SCALE'] = head_orig['SCALE']
					head_template['SCALE_E'] = head_orig['SCALE_E']
				except:
					pass

				try:
					head_template['H_OFFSET'] = head_orig['H_OFFSET']
					head_template['H_OFF_E'] = head_orig['H_OFF_E']
				except:
					pass

				# resample image 
				rf.regrid(files[i], template_map_name, resultimage = resamp_file, header = None)
			
			Final_names.append(resamp_file)
			print ()	

		print ()
	return Final_names

def generateGradientMaps_old (files, running = True):
	'''
	This function uses images to to get the gradient maps.
	
	Parameters
	----------
	file: list of strings
		Path and name of all the required files
	running: boolean
		If False, it will not create the new files but only 
		provide the name of the new files 

	Returns
	-------
	Final_names : List
		name of all the final files 

	'''

	print ()
	print ("GETING GRADIENT MAPS \nTotal number of files:", len(files))


	folder_name = 'gradient_maps/'
	print (len(folder))

	os.system('mkdir '+folder+folder_name)
	print ()

	Final_names = []

	for i in files:
		if running == True:
		
			# fix the header
			data, header = fits.getdata(i, header = True)

			try:
				header.remove(keyword='OFFSET')
				header.remove(keyword='OFFSET_E')
				header.remove(keyword='SCALE')
				header.remove(keyword='SCALE_E')
			except:
				print ()

			try:
				header.remove(keyword='H_OFFSET')
				header.remove(keyword='H_OFF_E')
				header.remove(keyword='PLN_CR1') 
				header.remove(keyword='PLN_CR2') 
			except:
				print()	
			
			header['EXTNAME'] = ('Gradient','name of this file')
			header['BUNIT']   = ('MJy/sr/pixel', 'Unit of the data') 

			
			# get the map
			gx, gy = np.gradient(data, 1,1)
			gradient= np.sqrt(gx**2+gy**2)
		

		file_name = i[0:len(folder)]+folder_name+i[len(folder):-5]+'.gradient.fits'
		
		if running == True:
			# save the map
			os.system('rm '+file_name)
			fits.writeto(file_name, gradient, header )
			current_res = beamsize[str(convolve_size)]
			final_res = np.sqrt(20**2+current_res**2)
			os.system('rm '+file_name[:-5]+'.conv.fits')
			rf.resconvolve(file_name, current_res, final_res )
			
		print ('File saved as:',file_name[:-5]+'.conv.fits')
	
		Final_names.append(file_name[:-5]+'.conv.fits')
	return Final_names

def generateGradientMaps (files, running = True):
	'''
	This function uses images to to get the gradient maps.
	
	Parameters
	----------
	file: list of strings
		Path and name of all the required files
	running: boolean
		If False, it will not create the new files but only 
		provide the name of the new files 

	Returns
	-------
	Final_names : List
		name of all the final files 

	'''

	print ()
	print ("GETING GRADIENT MAPS \nTotal number of files:", len(files))


	folder_name = 'gradient_maps/'
	print (len(folder))

	os.system('mkdir '+folder+folder_name)
	print ()

	Final_names = []

	for i in files:
		if running == True:
		
			# fix the header
			data, header = fits.getdata(i, header = True)

			try:
				header.remove(keyword='OFFSET')
				header.remove(keyword='OFFSET_E')
				header.remove(keyword='SCALE')
				header.remove(keyword='SCALE_E')
			except:
				print ()

			try:
				header.remove(keyword='H_OFFSET')
				header.remove(keyword='H_OFF_E')
				header.remove(keyword='PLN_CR1') 
				header.remove(keyword='PLN_CR2') 
			except:
				print()	
			
			header['EXTNAME'] = ('Gradient','name of this file')
			header['BUNIT']   = ('MJy/sr/pixel', 'Unit of the data') 

			
			# get the map

			#### ADD GRAD COMBINED HERE

			kernel = np.zeros([3,3])

			kernelxp = np.copy(kernel)
			kernelxp[1,0]=0
			kernelxp[1,1]=1
			kernelxp[1,2]=-1
			grad_xp = sig.convolve2d(data, kernelxp, mode='same')

			kernelxm = np.copy(kernel)
			kernelxm[1,0]=1
			kernelxm[1,1]=-1
			kernelxm[1,2]=0
			grad_xm = sig.convolve2d(data, kernelxm, mode='same')

			kernelyp = np.copy(kernel)
			kernelyp[0,1]=0
			kernelyp[1,1]=1
			kernelyp[2,1]=-1
			grad_yp = sig.convolve2d(data, kernelyp, mode='same')

			kernelym = np.copy(kernel)
			kernelym[0,1]=1
			kernelym[1,1]=-1
			kernelym[2,1]=0
			grad_ym = sig.convolve2d(data, kernelym, mode='same')

			a= np.sqrt(grad_xp**2 + grad_yp**2)
			b= np.sqrt(grad_xp**2 + grad_ym**2)
			c= np.sqrt(grad_xm**2 + grad_yp**2)
			d= np.sqrt(grad_xm**2 + grad_ym**2)

			gradient = ((a+b+c+d)/4)
		

		file_name = i[0:len(folder)]+folder_name+i[len(folder):-5]+'.gradient.fits'
		
		
		if running == True:
			# save the map
			os.system('rm '+file_name)
			fits.writeto(file_name, gradient, header )
			
		print ('File saved as:',file_name)
		Final_names.append(file_name)


	return Final_names	

def generateCPDMaps(wavelength, files, coverage_files, gradient_files, running =True):
	'''
	This function uses wavelengths and images to create the coverage predicted 
	dispersion maps.
	
	Parameters
	----------
	wavelengths : 1-D array
		array of all the wavelength in micron	
	file: list of strings
		Path and name of all the required intensity files
	coverage_files: list of strings
		Path and name of all the required coverage files
	gradient_files: list of strings
		Path and name of all the required gradient files
	running: boolean
		If False, it will not create the new files but only 
		provide the name of the new files 

	Returns
	-------
	Final_names : List
		name of all the final files 

	'''

	print ()
	print ("GETING UNCERTAINTY MAPS \nTotal number of files:", len(files))

	folder_name = 'uncertainty_maps/'

	os.system('mkdir '+folder+folder_name)	
	print()

	Final_names = []

	for num, i in enumerate(files):
		
		if running == True:
			# fix the header 
			header = fits.getheader(i)

			try:
				header.remove(keyword='OFFSET')
				header.remove(keyword='OFFSET_E')
				header.remove(keyword='SCALE')
				header.remove(keyword='SCALE_E')
			except:
				print ()

			try:
				header.remove(keyword='H_OFFSET')
				header.remove(keyword='H_OFF_E')
				header.remove(keyword='PLN_CR1') 
				header.remove(keyword='PLN_CR2') 
			except:
				print()	
			
			header['EXTNAME'] = ('Uncertainty','name of this file')
			header['BUNIT']   = ('MJy/sr', 'Unit of the data') 

			#get the map
			coverage = fits.getdata(coverage_files[num])
			gradient = fits.getdata(gradient_files[num])

			sbd_i = to_get_errors[str(wavelength[num])]['sbd_i']
			sbd_g = to_get_errors[str(wavelength[num])]['sbd_g']
			conv_i = to_get_errors[str(wavelength[num])]['conv_i']
			conv_g = to_get_errors[str(wavelength[num])]['conv_g']

			if num == 0: #for 160
				cpd = np.sqrt(((RIF_160**2 + WN_160**2))+(sbd_g*gradient*conv_g/np.sqrt(coverage))**2)
			else: # for all the other frequencies
				cpd = np.sqrt((sbd_i*conv_i/np.sqrt(coverage))**2+(sbd_g*gradient*conv_g/np.sqrt(coverage))**2)
		
		
		end_of_file_name = i[len(folder+field_name+'-'+str(wavelength[num])+'.offset.'):]
		
		file_name = i[0:len(folder)]+folder_name+field_name+'-'+str(wavelength[num])+'.cpd.'+end_of_file_name

		if running == True:
		
			# save the map
			os.system('rm '+file_name)
			fits.writeto(file_name, cpd, header )
		
		print ('File saved as:',file_name)		

		Final_names.append(file_name)


	return Final_names 

def generateMaskedMaps(result_folder):
	'''
	Masked the final maps for Common Science Region

	Parameters
	----------
	result_folder: strings
		Path to the final maps

	'''
	
	print ("Field:", field_name,"\n")

	path_to_mask = path_to_folder+'/grid'
	mask = '{0}/{1}_{2}.fits'.format(path_to_mask, field_name, regrid_size)

	path_to_map = result_folder
	maps=glob.glob('{0}/*.fits'.format(path_to_map))

	os.system('mkdir {0}/{1}_masked'.format(path_to_map,field_name))

	mask_data = fits.getdata(mask)
	mask_data[mask_data<=0.99] = float('nan')
	mask_data[mask_data>0.99] = 1

	for i in maps:
		data, head = fits.getdata(i, header=True)
		data_masked = data*mask_data
		fits.writeto(i[:-10]+'.fits', data_masked, head)
		os.system('mv {0} {1}/{2}_masked'.format(i[:-10]+'.fits',path_to_map,field_name))

	return


def create_radiance(temp, tau, beta_map, nu0):
	Temp = temp
	Tau = tau

	nu_0 = nu0#1*u.THz
	H,K, sig_SB = const.h , const.k_B, const.sigma_sb
	H,K, sig_SB = H.value , K.value , sig_SB.value

	beta_stuff = ((K*Temp/(H*nu_0))**beta_map)*((sp.zeta(4+beta_map)*sp.gamma(4+beta_map))/(sp.zeta(4)*sp.gamma(4)))
	Radiance = (Tau*sig_SB/np.pi)*(Temp**4)*beta_stuff
	return Radiance




	
