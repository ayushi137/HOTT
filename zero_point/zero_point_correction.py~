from __future__ import print_function
##############################################################################
'''
Herschel Offset using planck using the function script: required_functions
and offset_function

This script requires following functions and files:
initial_parameter.py
required_functions.py
colour_correction.py 
boundary.py --- Note: only need this to make the maps look more presentable
turbo_colormap.py 

To run:
run offset.py 'Planck field' 'Field' 'wavelength'
 
Example:
run offset.py TauLL TauS2 250 

By: Ayushi Singh
Last Modified: 16 March 2022
'''
##############################################################################
from astropy.io import fits
from astropy.wcs import WCS

import numpy as np 

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import rcParams
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

from sys import exit
import sys, os, datetime, time
import pandas as pd


# other scripts 
from initial_parameter import path_to_folder, beamsize
from initial_parameter import apply_cc, savecc, high_crop, percentage, gradient_crop, grad_percentage
from initial_parameter import low_crop, low_percentage, add_plane
from initial_parameter import residual_crop, sigma_clipping, edge_crop, res_lim
from initial_parameter import images, savefits, saveimage, write_offset
from initial_parameter import temperature, beta
import required_functions as rf
import boundary
import colour_correction as cc
from turbo_colormap import turbo_cmap


# alternative to tight layout
rcParams.update({'figure.autolayout': True})
plt.ion()

# Supress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=UserWarning) 

# wavelength, herschel and planck field name with file syntax
planck_name = sys.argv[1]
name = sys.argv[2]
wavelength = sys.argv[3]

extra = False # add extra stuff under the "if extra == True:" statement

print ('planck name: '+ planck_name)
print ('field name: '+ name)
print ('wavelength: '+ wavelength + ' micron')
print ()
print ('adding plane?:', add_plane)
print ()


#####################################################################################################
################################## VARIOUS FILE NAMES ###############################################
#####################################################################################################

main_folder = path_to_folder

# input file
herschel_file = main_folder+'/herschel_data/{0}/{0}-{1}.image.fits'.format(name, wavelength)
planck_file = main_folder+'/Models/{0}_HFI_model_{1}.fits'.format(planck_name, wavelength)

# noise maps
noise_file_orig = main_folder+'/herschel_data/{0}/sigma_maps/{0}-{1}.sigma.fits'.format(name, wavelength)
noise_file_temp = main_folder+'/herschel_data/{0}/sigma_maps/{0}-250.sigma.fits'.format(name, wavelength)

# convolved and resampled file
conv_file = main_folder+'/offset_maps/{0}-{1}.image.conv.fits'.format(name, wavelength)
resamp_file = main_folder+'/offset_maps/{0}-{1}.image.conv.resamp.fits'.format(name, wavelength)
her_corr = main_folder+'/offset_maps/{0}-{1}'.format(name, wavelength)
masked_mask_name = main_folder+'/offset_mask/{0}-{1}.masked_mask.fits'.format(name, wavelength)

# output image and herschel corrected file
dateToday = (datetime.date.today()) 
os.system('mkdir {0}/offset_images/offsets-{1}'.format(main_folder,str(dateToday)))
os.system('mkdir {2}/herschel_data/{0}/planes'.format(name, wavelength,main_folder))

output_combine = main_folder+'/offset_images/offsets-{2}/offsets-{0}-{1}.png'.format(name, wavelength, str(dateToday)) 
output_file = main_folder+'/herschel_data/{0}/{0}-{1}.offset.fits'.format(name, wavelength)
output_plane = main_folder+'/herschel_data/{0}/planes/{0}-{1}.plane.fits'.format(name, wavelength)


# masking file for the plnack grid
mask_file = main_folder+'/grid/{0}_planck.fits'.format(name)

# planck dustmodel map
temp_map = main_folder+'/Planck/{0}{1}'.format(planck_name, temperature)
beta_map = main_folder+'/Planck/{0}{1}'.format(planck_name, beta)

# saving cc map
cc_map = main_folder+'/cc_maps/cc-{0}-{1}.fits'.format(name, wavelength) 

# saving the offset values in a txt file
offset_csv = main_folder+'/offset_values/offset_values_'+str(dateToday)+'.csv'
#####################################################################################################
######################################## FUNCTIONS ##################################################
#####################################################################################################
def applying_cc(intensity_map, wavelength,temp_map, beta_map):

	"""
	Apply colour correction to intensity at a given wavelength using the temperature and beta values. 
	It uses colour_correction.py

	Parameters
	----------
	intensity_map : 2D image
		intensity map
	wavelength : int or float
		wavelength value
	temp_map : string
		path to temperature map
	beta_map : string
		path to beta map		

	Returns
	-------
	intensity_map : 2D image
		colour corrected intensity map
	colour_corrections : 2D image
		colour correction image

	"""	
	# get the temperature and beta from intensity_map maps
	temperature = fits.getdata(temp_map)
	beta = fits.getdata(beta_map)

	# get colour correction value from the colour_correction.py
	colour_corrections = np.zeros([np.shape(temperature)[0],np.shape(temperature)[1]])
	for i in range(np.shape(beta)[0]):
		for j in range(np.shape(beta)[1]):
			colour_corrections[i][j] = cc.colourcorrection(temperature[i][j], beta[i][j], wavelength = int(wavelength))
	        
	# correct intensity_map 
	intensity_map = intensity_map/colour_corrections

	return intensity_map, colour_corrections

def cropping(name, data):
	"""
    This crops the images if proper values are given to file boundary.py

    Parameters
    ----------
    name : sting
    	name given to the region in this code 
    data : 2D image
    	data map	

    Returns
    -------
    data : 2D image
    	cropped data map
   	"""	

	top,left,bottom,right = boundary.boundaries(name, 'grid')						

	data = data[bottom:top,left:right]

	print ('Done cropping.\n')

	return data
		
def leastfit(x,y, write):
	"""
    Use chi reduced method to get linear line fit

    Parameters
    ----------
    x : float
    	the independent points
    y : float
    	the dependent points
    write : boolean
    	if True, then it prints the result with error
    	if False, then it doesn't print anything	

    Returns
    -------
    polyfit : 1-D array 
        slope and the offset value
    resudial: 1-D array     
    	resudial between he line and the data points 
   	"""	

	ma =np.array([[np.sum(x**2),np.sum(x)],[np.sum(x),len(x)]])
	mc =np.array([[np.sum(x*y)],[np.sum(y)]])

	mai = np.linalg.inv(ma)
	md = np.dot(mai,mc)

	mfit = md[0,0]
	cfit = md[1,0]
	
	n = len(x)
	sum_x = np.sum(x)
	sum_x2 = np.sum(x**2)

	variance = (1.0/(len(x)-2))*np.sum(y-x*mfit-cfit)**2.0
	err_m = np.sqrt(n*variance/((n*sum_x2)-(sum_x)**2))
	err_c = np.sqrt(variance*sum_x2/((n*sum_x2)-(sum_x)**2))

	residual = y - (mfit*x+cfit)

	polyfit = [mfit , cfit]
	if write == True:	
		print ()
		print ('mfit:',mfit)
		print ('mfit error:',err_m)
		print ('cfit:', cfit)
		print ('cfir error:', err_c)
		print ('variance:', variance)
		print ()
	
	return polyfit 

def peval(x, p):  	
	return ((p[0]*x)+p[1])

def residuals (p,y,x, peval):
	return (y) - peval(x,p)	

def leastsqfit(p0,x,y):
	"""
    Use chi reduced method to get linear line fit using the python module

    Parameters
    ----------
    p0 : list
    	list of the parameters
    x : float
    	the independent points
    y : float
    	the dependent points	

    Returns
    -------
    polyfit : truple
        solution, error on the solution, status, covarianve matrix resulted from
        the fitting program 
   	"""	

	polyfit = leastsq(residuals,p0,args=(y,x, peval), full_output= True)

	y_final = peval(x,polyfit[0])
	chi2 = np.sum((y - y_final)**2/ y_final)# ((dy)**2))
	resi = (residuals(polyfit[0],y,x,peval))
	dof = len(y)-len(p0)
	chi_re2 = chi2/dof # residual variance
	cov = polyfit[1] * chi_re2
	cov_xy = cov[0][1]
	cov_x = np.sqrt(cov[0][0])
	cov_y = np.sqrt(cov[1][1])
	r =cov_xy/cov_x*cov_y

	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ("The inital parameter (p0) we used is:\n" , p0)
	print ("What we get as a parameter:\n" , polyfit[0])

	if polyfit[4] == 1 or polyfit[4] == 2 or polyfit[4] == 3 or polyfit[4] == 4: 
	    print ("It converges.")
	else:
	    print ("It does not converge.")

	print ("The Chi square is: \t\t\t" , round(chi2,2))
	print ("The Chi-reduced square is: \t\t" , round(chi_re2,2))
	print ()
	print ("Cov_xy:" , round(cov_xy,4) , "\nCov_x: " , round(cov_x,4) , "\nCov_y: " , round(cov_y,4))
	print ("error_x:" , round(np.sqrt(cov_x),4) , "\nerror_y:" , round(np.sqrt(cov_y),4))
	print ("Sample coefficient of linear correlation (r): " , round(pearsonr(x,y)[0],2))
	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ()
	return polyfit, cov_x, cov_y


def f_min(X,p):
    plane_xyz = p[0:3]
    if add_plane == True:
    	distance = p[2]*(p[0]*X.T[0] + p[1]*X.T[1] +X.T[2])+ p[3]
    else:
    	distance = p[0]*(X.T[2])+ p[1]	
    return distance 
    
def residuals( params, Y_data, X, f_min):
    return Y_data - f_min(X, params)

def fitting(y, x, p0, X_mid, Y_mid, plotting = True): 

	"""
    Use chi reduced method to get fit for a plane, scale and zero-point correction.
    We first fit y on x then we change it to x on y. 

    Parameters
    ----------
    x : float
    	the independent points
    y : float
    	the dependent points	
    p0 : list
    	list of the parameters
    X_mid : int
    	middle of the data in X
    Y_mid : int
    	middle of the data in Y		

    Returns
    -------
    solution : array 
    	resulting parameters
    plane : 2D aray
    	resulting plane
    error : array
    	errors on the parameter
   	"""	  

	print ('********* We fit Herschel on Planck *********\n')  
	XYZ = y   ### X axis [Herchel]
	Y_all = x   #### y axis [Planck]
	#print (np.shape(XYZ))
	Z = XYZ.ravel()

	mask = (np.isnan(Z) == False)
	Z = Z[mask]
	X = np.arange(np.shape(XYZ)[1])
	Y = np.arange(np.shape(XYZ)[0])
	XX,YY= np.meshgrid(X,Y)

	XX = XX - X_mid
	YY = YY - Y_mid

	#XX = XX - (np.shape(XYZ)[1]//2)
	#YY = YY - (np.shape(XYZ)[0]//2)

	xs, ys = XX.ravel()[mask], YY.ravel()[mask]
	zs = XYZ.ravel()[mask]

	Y_data = Y_all.ravel()[mask]

	XYZ = (np.c_[xs,ys,zs])

	#p0 = [0.5, 0.5, polyfit[0], polyfit[1]] 

	final = leastsq(residuals, p0, args=(Y_data, XYZ, f_min), full_output= True)

	if add_plane == True:

		y_final = f_min(XYZ,final[0])
		chi2 = np.sum((Y_data - y_final)**2/ y_final)# ((dy)**2))
		dof = len(Y_data)-len(p0)
		chi_re2 = chi2/dof # residual variance
		cov = final[1] * chi_re2
		cov_xy = cov[0][1]
		cov_ab = cov[2][3]
		error_x = np.sqrt(cov[0][0])
		error_y = np.sqrt(cov[1][1])
		error_a = np.sqrt(cov[2][2])
		error_b = np.sqrt(cov[3][3])
		r =cov_ab/error_a*error_b

		sol = final[0]
		print (sol)

	else:

		y_final = f_min(XYZ,final[0])
		chi2 = np.sum((Y_data - y_final)**2/ y_final)# ((dy)**2))
		dof = len(Y_data)-len(p0)
		chi_re2 = chi2/dof # residual variance
		cov = final[1] * chi_re2
		cov_xy = cov[0][1]
		#cov_ab = cov[2][3]
		error_x = np.sqrt(cov[0][0])
		error_y = np.sqrt(cov[1][1])
		error_a = 0#np.sqrt(cov[2][2])
		error_b = 0#np.sqrt(cov[3][3])
		#r =cov_ab/error_a*error_b

		sol = np.array([0, 0, final[0][0], final[0][1]]) 

	print (sol)

	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ("The inital parameter (p0) we used is:\n" , p0)
	print ("What we get as a parameter:\n" , sol)

	if final[4] == 1 or final[4] == 2 or final[4] == 3 or final[4] == 4: 
	    print ("It converges.")
	else:
	    print ("It does not converge.")

	print ("The Chi square is: \t\t\t" , round(chi2,2))
	print ("The Chi-reduced square is: \t\t" , round(chi_re2,2))
	print ()
	print ("error_x:" , round(error_x,4) , "\nerror_y:" , round(error_y,4))
	print ("error_a:" , round(error_a,4) , "\nerror_b:" , round(error_b,4))
	print ("Sample coefficient of linear correlation (r): " , round(pearsonr(Y_data, XYZ[:,2])[0],2))
	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ()


	print ('********* We change it to Planck on Herschel  *********\n')

	solution = np.copy(sol)
	solution[2] = 1/sol[2]
	solution[3]=   sol[3]/sol[2]

	error_x = error_x/sol[2]
	error_y = error_y/sol[2]

	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ("What we get as a parameter:\n" , solution)
	print ()
	print ("error_x:" , round(error_x,4) , "\nerror_y:" , round(error_y,4))
	print ("error_a:" , round(error_a,4) , "\nerror_b:" , round(error_b,4))
	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ()

	xx, yy = XX,YY
	z_plane = (solution[0]*xx+solution[1]*yy)

	if plotting == True:
		from mpl_toolkits.mplot3d import Axes3D
		fig = plt.figure()
		ax = fig.gca(projection='3d')

		ax.view_init(elev=6, azim=-11)
		ax.scatter(XYZ[:, 0], XYZ[:, 1], Y_data -(XYZ[:, 2]+solution[3])/solution[2], c= Y_data -(XYZ[:, 2]+solution[3])/solution[2], cmap='inferno')
		ax.plot_surface(xx, yy, z_plane, alpha=0.5, cmap = turbo_cmap)
		plt.title('{0}, {1}, {2}, {3}'.format(round(solution[0],4),round(solution[1],4),round(solution[2],4),round(solution[3],4)))

	errors = np.array([error_x, error_y, error_a, error_b])
	return solution, z_plane, errors

def high_cropping_fun(cuttoff_map, percentage):
	mask_pixel= (np.isnan(cuttoff_map)==False)
	no_nan = cuttoff_map[mask_pixel]

	# make histogram
	bin_num = 1000
	
	k = 0 

	while k <= 0:
		hist, bin_val = np.histogram(no_nan, bin_num)
		
		i = 0
		j = 0
		
		all_hist = hist[::-1]
		all_bin = bin_val[::-1]
		while i < np.sum(hist)*(percentage/100.):
			i += all_hist[j]
			j += 1
	
	
		k = bin_num - j 
		#print (j, k)

		bin_num = bin_num+100
			
	return bin_val, k

def low_cropping_fun(low_cuttoff_map, low_percentage):
	mask_pixel= (np.isnan(low_cuttoff_map)==False)
	no_nan = low_cuttoff_map[mask_pixel]

	# make histogram
	bin_num = 1000
	
	k = 0 

	while k <= 0:
		hist, bin_val = np.histogram(no_nan, bin_num)
		
		i = 0
		j = 0
		
		all_hist = hist
		all_bin = bin_val
		while i < np.sum(hist)*(low_percentage/100.):
			i += all_hist[j]
			j += 1
	
	
		k = j 
		#print (j, k)

		bin_num = bin_num+100
	return bin_val, k 

# Equation for Gaussian
def f(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))
#####################################################################################################
################################## MAIN PROGRAM #####################################################
#####################################################################################################
# start timer 
start_timer = time.time()

##### Get the mask file for the planck grid
mask = fits.getdata(mask_file)
mask[mask <= 0.01] = float('nan')
mask[mask > 0.01] = 1
#mask = mask[0:211,141:423]

##### convolve and regrid the herschel image to planck resolution and pixelsize 
# remove if these files are already there
os.system("rm {0}".format(conv_file))
os.system("rm {0}".format(resamp_file))
# convolve the image to planck resolution of 5'
rf.resconvolve(herschel_file,beamsize[wavelength], beamsize['planck'] ,resultimage = conv_file)
# resample image to planck pixel size
rf.reproject(conv_file, planck_file, resultimage = resamp_file, header = None)

##### get the proper planck and herchel resized image
herschel, herschelhead= fits.getdata(resamp_file, header = True)
planck, planckhead = fits.getdata(planck_file, header = True)
print ('Herschel and Planck images have been loaded.\n')

herschel_orig = np.copy(herschel)

##### mask them to get the only common science region
planck = planck/mask
herschel = herschel/mask

herschel_orig = np.copy(herschel)

## herschel data in the native grid ####
herschel_native ,head = fits.getdata(herschel_file, header = True)


######################### applying colour correction to planck model #########################

if apply_cc == True:
	print ('Start the colour correction script.')
	planck, colour_corrections = applying_cc(planck, wavelength, temp_map, beta_map)	
	print ('Finish colour correcting Planck data.\n')


mean_before_crop = np.nanmean(planck)
median_before_crop = np.median(planck[~np.isnan(planck)])
median_before_crop_her= np.median(herschel[~np.isnan(herschel)])


######################### remove high values from herschel that won't agree well with planck ####################

if high_crop == True:
	# remove some percent value values 
	print ('###########################################################################')
	print ('Start removing top ', percentage ,r'% of the data.')

	cuttoff_map = np.copy(planck)

	bin_val, k = high_cropping_fun(cuttoff_map, percentage)
	
	planck[cuttoff_map > bin_val[k]] = float('nan')
	herschel[cuttoff_map > bin_val[k]] = float('nan')

	print ('Finish removing the higher points.\n')
	print ('###########################################################################')

######################### remove low  values from herschel that won't agree well with planck ##################
if low_crop == True:
	# remove some percent value values 
	print ('###########################################################################')
	print ('Start getting bottom ', low_percentage ,r'% of the planck data.')

	low_cuttoff_map = np.copy(planck)

	bin_val, k = low_cropping_fun(low_cuttoff_map, low_percentage)
	
	low_mean = np.nanmean(planck[low_cuttoff_map > bin_val[j]])
	planck[low_cuttoff_map < bin_val[j]] = float('nan')
	herschel[low_cuttoff_map < bin_val[j]] = float('nan')

	print ('Finish removing the lower points.\n')
	print ('###########################################################################')


####### gradient cutoff

if gradient_crop == True:

	print ('###########################################################################')
	print ('Start removing top ', grad_percentage ,r'% of the gradient.')
	gradx,grady = np.gradient(herschel_orig, 1,1)
	grad = np.sqrt(gradx**2+grady**2)

	cuttoff_map = np.copy(grad)

	bin_val, k = high_cropping_fun(cuttoff_map, grad_percentage)
	
	planck[cuttoff_map > bin_val[k]] = float('nan')
	herschel[cuttoff_map > bin_val[k]] = float('nan')

	print ('Finish removing the higher gradient points.\n')
	print ('###########################################################################')

########################## remove low  values from herschel that won't agree well with planck ##################

if residual_crop == True:
	print ('###########################################################################')
	print ('Start removing very high and low resudial points by sigma clipping of',sigma_clipping,'sigma.')
	mask_pixel2= (np.isnan(herschel)==False) 

	x_meas1 = planck[mask_pixel2]
	y_meas1 = herschel[mask_pixel2]

	maxp = np.max(x_meas1)

	top,left,bottom,right = boundary.boundaries(name, 'grid')

	if top == None:
		top = np.shape(herschel)[0]-1

	if bottom == None:
		bottom = 0

	if left == None:
		left = 0

	if right == None:
		right = np.shape(herschel)[1]-1

	X_half = (right - left)//2
	Y_half = (top-bottom)//2

	X_mid = left + X_half
	Y_mid = bottom + Y_half

	plt.figure(figsize = (9,8))

	plt.subplot(221)
	if edge_crop == True:
		planck_hold = cropping(name, planck)
		herschel_hold = cropping(name, herschel)

	polyfit = leastfit(x_meas1,y_meas1, write = True)
	residual_z = planck_hold-((herschel_hold-polyfit[1])/polyfit[0])
	residuals_z = residual_z[np.isnan(residual_z)==False]

	planckhead['OFFSET'] = (polyfit[1], '[MJy/sr] zero-point, from Planck dust models')
	planckhead['SCALE'] = (polyfit[0], 'scale factor, drived from Planck dust models')
	
	os.system("rm {0}.zeropoint.scale.fits".format(her_corr))
	fits.writeto(her_corr+'.zeropoint.scale.fits', ((herschel_orig-polyfit[1])/polyfit[0]), planckhead)

	os.system("rm {0}.zeropoint.fits".format(her_corr))
	fits.writeto(her_corr+'.zeropoint.fits', ((herschel_orig-polyfit[1])), planckhead)

	totalpixel = np.count_nonzero(~np.isnan(herschel)) 
	bin_num = int(totalpixel/10)
	n1, bins1, patches = plt.hist(residuals_z.ravel(), bins=bin_num, color='darkorange',density=True, alpha = 0.7,align='right') 
	popt1, pcov1 = curve_fit(f, bins1[1:], n1)
	y1 = f(bins1[1:], *popt1)
	plt.plot(bins1[1:], y1, c='k', linewidth=1, label = '{2}: $\mu = ${0}, $\sigma = ${1}'.format(round(popt1[1],4), round(popt1[2],4), wavelength))
	plt.legend()
	plt.title('Without Plane')

	plt.subplot(222)
	plt.imshow(residual_z, origin = 'lower', vmin = maxp*res_lim*-1, vmax = maxp*res_lim, cmap = turbo_cmap); plt.colorbar()
	plt.title('{0} at Wavelength {1}'.format(name, wavelength) + r'$\rm{\mu m}$')


	if add_plane == True:
		p0=[0.5, 0.5, polyfit[0], polyfit[1]]

	else:
		p0=[polyfit[0], polyfit[1]]

	solution, z_plane, errors = fitting(herschel,planck,p0, X_mid, Y_mid, plotting = False)

	if edge_crop == True:
		planck_hold = cropping(name, planck)
		herschel_hold = cropping(name, herschel)
		z_plane_hold = cropping(name, z_plane)

	resi = (planck-((herschel+z_plane+solution[3])/solution[2]))

	planckhead['OFFSET'] = (solution[3], '[MJy/sr] zero-point, from Planck dust models')
	planckhead['PLANEX'] = (solution[0], 'scale factor in x, to correct the plane')
	planckhead['PLANEY'] = (solution[1], 'scale factor in x, to correct the plane')
	planckhead['SCALE'] = (solution[2], 'scale factor, drived from Planck dust models')

	os.system("rm {0}.zeropoint.plane.scale.fits".format(her_corr))
	fits.writeto(her_corr+'.zeropoint.plane.scale.fits', ((herschel_orig+z_plane+solution[3])/solution[2]), planckhead)

	os.system("rm {0}.zeropoint.plane.fits".format(her_corr))
	fits.writeto(her_corr+'.zeropoint.plane.fits', ((herschel_orig+z_plane+solution[3])), planckhead)
	
	plt.subplot(223)
	residual_z = (planck_hold-((herschel_hold+z_plane_hold+solution[3])/solution[2]))
	residuals_z = residual_z[np.isnan(residual_z)==False]
	n1, bins1, patches = plt.hist(residuals_z.ravel(), bins=bin_num, color='teal',density=True, alpha = 0.7,align='right') 
	popt1, pcov1 = curve_fit(f, bins1[1:], n1)
	y1 = f(bins1[1:], *popt1)
	plt.plot(bins1[1:], y1, c='k', linewidth=1, label = '{2}: $\mu = ${0}, $\sigma = ${1}'.format(round(popt1[1],4), round(popt1[2],4), wavelength))
	plt.legend()
	plt.title('With Plane')

	plt.subplot(224)
	plt.imshow(residual_z, origin = 'lower', vmin = maxp*res_lim*-1, vmax = maxp*res_lim, cmap = turbo_cmap); plt.colorbar()

	residual_cutoff = abs(popt1[2])*sigma_clipping
	print ('Residual masking is from +/-', residual_cutoff)
	
	resi[resi > residual_cutoff] = float('nan')
	resi[resi < -1*residual_cutoff] = float('nan')
	resi[resi <= residual_cutoff] = 1.0
	resi[resi >= -1*residual_cutoff] = 1.0

	planck = planck/resi
	herschel = herschel/resi

	if saveimage == True:	
		print ('\nSaving png file.')
		plt.savefig('{0}'.format(output_combine[:-4]+'.before.png'))
	
	print ('Finish removing very high and low resudial points. \n')
	print ('###########################################################################')

masked_mask = np.copy(herschel)
masked_mask[np.isnan(herschel)==False] = 1
masked_mask[np.isnan(herschel)==True] = 0

os.system('rm {0}'.format(masked_mask_name))
fits.writeto(masked_mask_name, masked_mask, planckhead )
# final cropping for make it look readable
#if edge_crop == True:
#	planck = cropping(name, planck)
#	herschel = cropping(name, herschel)
	#z_plane = cropping(name, z_plane)
######################### Fitting #############################################
# NOTE: we are correlating y on x that is herschel on planck

# remove all nan
print ('###########################################################################')
mask_pixel3= (np.isnan(herschel)==False) & (np.isnan(planck)==False)

x_meas1 = planck[mask_pixel3]
y_meas1 = herschel[mask_pixel3]

# this is used to create the line of fit
start = np.min([0,np.min(x_meas1)])
test = np.arange(0,np.max(x_meas1),0.01)

# fitting the line and getting the slope and the offset	estimate
polyfit = leastfit(y_meas1, x_meas1, write = True)
#polyfit = np.polyfit(y_meas1,x_meas1, 1)

top,left,bottom,right = boundary.boundaries(name, 'grid')

if top == None:
	top = np.shape(herschel)[0]-1

if bottom == None:
	bottom = 0

if left == None:
	left = 0

if right == None:
	right = np.shape(herschel)[1]-1

X_half = (right - left)//2
Y_half = (top-bottom)//2

X_mid = left + X_half
Y_mid = bottom + Y_half

p0=[0.5, 0.5, polyfit[0], polyfit[1]]

solution, z_plane, errors = fitting(herschel,planck,p0, X_mid, Y_mid)

error_x, error_y, error_a, error_b = errors

# plane to save using the whole planck image 
plane_save = z_plane#(solution[0]*xx+solution[1]*yy)
print ('\nSaving fits file for plane and reggriding it.')

planckhead['OFFSET'] = (solution[3], '[MJy/sr] zero-point, from Planck dust models')
planckhead['PLANEX'] = (solution[0], 'scale factor in x, to correct the plane')
planckhead['PLANEY'] = (solution[1], 'scale factor in x, to correct the plane')
planckhead['SCALE'] = (solution[2], 'scale factor, drived from Planck dust models')

os.system("rm {0}".format(output_plane))
fits.writeto('{0}'.format(output_plane), plane_save, planckhead)

os.system("rm {0}.offset.scale.fits".format(her_corr))
fits.writeto(her_corr+'.offset.scale.fits', ((herschel_orig+z_plane+solution[3])/solution[2]), planckhead)

os.system("rm {0}.offset.fits".format(her_corr))
fits.writeto(her_corr+'.offset.fits', ((herschel_orig+z_plane+solution[3])), planckhead)


# to make a plane at herschel pixel scale 
wp = WCS(planckhead)
wh = WCS(head)
planck_scale = planckhead['CDELT2']
herschel_scale = head['CDELT2']

scale_x = solution[0]*herschel_scale/planck_scale
scale_y = solution[1]*herschel_scale/planck_scale

error_x_new = error_x*herschel_scale/planck_scale
error_y_new = error_y*herschel_scale/planck_scale

ra_mid, dec_mid = wp.wcs_pix2world(X_mid, Y_mid, 0)

X_mid_h, Y_mid_h = wh.wcs_world2pix(ra_mid, dec_mid, 0)

X_mid_h, Y_mid_h = float(X_mid_h), float(Y_mid_h)

print ('Plane centre pixel in Planck:\t', X_mid, Y_mid)
print ('Plane centre coordinate:\t', ra_mid, dec_mid)
print ('Plane centre pixel in Native:\t', X_mid_h, Y_mid_h)
print ()
print ('Pixel scale in Planck in arcsec:', planck_scale*3600)
print ('Pixel scale in Native in arcsec:', herschel_scale*3600)
print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print ()

map_grid = np.copy(herschel_native)
Xh = np.arange(np.shape(map_grid)[1])
Yh = np.arange(np.shape(map_grid)[0])
XXh,YYh= np.meshgrid(Xh,Yh)

XXh = XXh - X_mid_h#- (np.shape(map_grid)[1]//2)
YYh = YYh - Y_mid_h#- (np.shape(map_grid)[0]//2)

xx, yy = XXh,YYh
plane = (scale_x*xx+scale_y*yy) # plane at herschel pixel scale

print ('The parameters are:\t', round(scale_x,6), round(scale_y,6), round(solution[2],6), round(solution[3],6))
print ('The parameter errors:\t', round(error_x_new,6), round(error_y_new,6), round(error_a,6), round(error_b,6))

print ('###########################################################################')
#plane_regrid = rf.regrid(output_plane, herschel_file, savefile = False)

#plane = plane_regrid#-np.nanmean(plane_regrid)

#os.system("rm {0}".format(output_plane[:-5]+'.resamp.fits'))
#fits.writeto('{0}'.format(output_plane[:-5]+'.resamp.fits'), plane_regrid, head)

os.system("rm {0}".format(output_plane[:-5]+'.native.fits'))
fits.writeto('{0}'.format(output_plane[:-5]+'.native.fits'), plane, head)

# final cropping for make it look readable
if edge_crop == True:
	planck = cropping(name, planck)
	herschel = cropping(name, herschel)
	z_plane = cropping(name, z_plane)

if extra == True:
	herschel = herschel[:,:60]
	planck = planck[:,:60]

######################### Plotting #############################################
# this is all plotting from here on
try:
		
	# MAKING THE IMAGE

	# creating the figure size
	plt.figure(figsize=(15,7))
	
	# SCATTER PLOT

	mask_pixel4 = (np.isnan(herschel)==False) & (np.isnan(planck)==False)

	x_meas1 = planck[mask_pixel4]
	y_meas1 = (herschel+z_plane)[mask_pixel4] ###### DOUBBLE CHECK THE SIGN
	test = np.arange(0,np.max(x_meas1),0.01)
	cpad = 0.01

	# Setting up
	#a = np.round(polyfit1[0],2)
	#b = np.round(polyfit1[1],2)

	a = solution[2]
	b = solution[3]	
	#labelfit = '{2} \n= {0} Planck + {1}'.format(a,b, name)
	labelfit = 'Herschel \n= {0} Planck - {1}'.format(np.round(a,2),np.round(b,2), name)
	x_label = 'Planck-based Dust Model [MJy/sr]'
	y_label = 'Herschel Intensity incl. plane [MJy/sr]'.format(name)
	title = '{0} at {1}'.format(name, wavelength) + r'$\rm{\mu m}$'

	def peval2(x, p):  	
		return ((p[0]*x)-p[1])

	# Plotting 
	if images == True:
		plt.subplot2grid((2,3),(0, 0),colspan=1, rowspan=2)	

	H, xedges, yedges = np.histogram2d(x_meas1,y_meas1, bins=70)
	H = np.rot90(H)
	H = np.flipud(H)	
	Hmask = np.ma.masked_where(H==0, H)
	plt.pcolormesh(xedges,yedges, Hmask, vmin = 0, vmax = 30, cmap = turbo_cmap)
	cbar=plt.colorbar(use_gridspec=True, pad = cpad)
	cbar.set_label('Pixel density',rotation=270,labelpad=20,fontsize=15)
	cbar.ax.tick_params(labelsize=15)

	plt.plot(test,(peval2(test,[a,b])),color='w',markevery=10,linewidth = 5, alpha = 0.5)
	plt.plot(test,(peval2(test,[a,b])),color='maroon',markevery=10,linewidth = 3,label=labelfit)
	#plt.plot(test,(peval2(test,[1.0,(b/a)])), color='b', alpha=0.5, linewidth = 1 )
	#plt.plot([0,np.max(x_meas1)],[-b/a, (peval2(np.max(x_meas1),[a,b]))], color='k', alpha=0.5, linewidth = 1 )
	#plt.scatter(0, -b/a, marker = 'x', color = 'k', linewidths = 3, label= '{0}'.format(round(-b/a,2)))

	#plt.plot(test,(f_min(test,sol)),color='maroon',markevery=10,linewidth = 1,label=labelfit)
	#plt.legend(loc=2,ncol=1,prop={'size':9})
	plt.legend(loc=2,ncol=1,prop={'size':15}, handlelength=1)
	plt.title(title,fontsize=15)
	plt.xlabel(x_label,fontsize=15)
	plt.ylabel(y_label,fontsize=15)
	plt.tick_params(labelsize=15)

	
	if images == True:
		# IMAGES OF TWO MAPS AND RESIDUAL

		line = r'${\mid}$'
		maxp = np.max(x_meas1)
		minp = np.min(x_meas1)
		low_lim = minp #- (minp%10)
		high_lim = maxp #- (maxp%10) #+ 10
		
		
		# vmin and vmax are lower and upper limit of the colorbar. For planck it was defined earlier 
		# and for hersechel it was calualted using the offset and scale 

		# PLANCK
		plt.subplot2grid((2,3),(0, 1))
		plt.imshow(planck,origin='lower',vmin =low_lim, vmax=high_lim, cmap = turbo_cmap)
		plt.title('Planck-based Dust Model '.format(title, wavelength),fontsize=15)
		plt.tick_params(labelsize=15)
		#plt.xlabel('Pixel [x]')
		#plt.ylabel('Pixel [y]')
		#divider = make_axes_locatable(plt.gca())
		#cax = divider.append_axes("right"+ "5%", pad="3%")
		cbar=plt.colorbar(use_gridspec=True, pad = cpad)
		cbar.ax.tick_params(labelsize=15)
		#cbar=plt.colorbar()#(cax=cax)
		cbar.set_label('[MJy/sr]',rotation=270,labelpad=15,fontsize=15)

		
		# HERSCHEL MAP
		plt.subplot2grid((2,3),(0, 2))
		#plt.imshow(herschel,origin='lower',vmin =(a*low_lim+b), vmax=(a*high_lim+b), cmap = turbo_cmap)
		plt.imshow((herschel+z_plane+b)/a,origin='lower',vmin =(low_lim), vmax=(high_lim), cmap = turbo_cmap)
		plt.title('Planck-corrected Herschel'.format(name),fontsize=15)
		plt.tick_params(labelsize=15)
		#plt.xlabel('Pixel [x]')
		#divider = make_axes_locatable(plt.gca())
		#cax = divider.append_axes("right", "5%", pad="3%")
		cbar=plt.colorbar(use_gridspec=True, pad = cpad)
		cbar.ax.tick_params(labelsize=15)
		#cbar=plt.colorbar()#(cax=cax)
		cbar.set_label('[MJy/sr]',rotation=270,labelpad=15,fontsize=15)

		# RESIDUAL 

		# Residuals are calculated using (Planck - (Herschel to Planck))
		#title = 'Subtracted Residual'
		title = 'Fractional Residual'
		plt.subplot2grid((2,3),(1, 1))
		#plt.imshow((planck-(herschel-b)/a)/planck,origin='lower',vmin = res_lim*-1, vmax=res_lim, cmap = turbo_cmap)
		#vmin = high_lim*res_lim*-1, vmax=high_lim*res_lim 
		plt.imshow((planck-(herschel+z_plane+b)/a)/planck,origin='lower',vmin = res_lim*-1, vmax=res_lim,  cmap = turbo_cmap)
		plt.title(r'{0}'.format(title, wavelength),fontsize=15)
		plt.tick_params(labelsize=15)
		#plt.xlabel('Pixel [x]')
		#plt.ylabel('Pixel [y]')
		#divider = make_axes_locatable(plt.gca())
		#cax = divider.append_axes("right", "5%", pad="3%")
		cbar=plt.colorbar(use_gridspec=True, pad = cpad)
		cbar.ax.tick_params(labelsize=15)
		#cbar=plt.colorbar()#(cax=cax)
		cbar.set_label('[Fraction]',rotation=270,labelpad=15,fontsize=15)	

		# Residuals are calculated using (Planck - (Herschel to Planck + plane))
		title = 'Residual '
		plt.subplot2grid((2,3),(1, 2))
		#plt.imshow((planck-(herschel-b)/a),origin='lower',vmin = high_lim*res_lim*-1, vmax=high_lim*res_lim , cmap = turbo_cmap)
		plt.imshow((planck-(herschel+b+z_plane)/a),origin='lower',vmin = high_lim*res_lim*-1, vmax=high_lim*res_lim ,  cmap = turbo_cmap)
		plt.title(r'{0}'.format(title, wavelength),fontsize=15)
		plt.tick_params(labelsize=15)
		#plt.xlabel('Pixel [x]')
		#divider = make_axes_locatable(plt.gca())
		#cax = divider.append_axes("right", "5%", pad="3%")
		cbar=plt.colorbar(use_gridspec=True, pad = cpad)
		cbar.ax.tick_params(labelsize=15)
		#cbar=plt.colorbar()#(cax=cax)
		cbar.set_label('[MJy/sr]',rotation=270,labelpad=15,fontsize=15)	

	
	#herschel_native = (herschel_native+((plane+b)/a))
	#print ('WITH THE SCALE')
	herschel_native = (herschel_native+((plane+b)))
	print ('WITHOUT THE SCALE')


	head['OFFSET'] = (b, '[MJy/sr] zero-point, from Planck dust models')
	head['OFFSET_E'] = (error_b, '[MJy/sr] zero-point error, error on zero-point')
	head['SCALE'] = (a, 'scale factor, drived from Planck dust models')
	head['SCALE_E'] = (error_a, 'scale factor error, error on scale factor value')
	head['PLANEX'] = (scale_x, 'scale factor in x, to correct the plane')
	head['PLANEX_E'] = (error_x_new, 'scale factor error, error on the scale factor')
	head['PLANEY'] = (scale_y, 'scale factor in x, to correct the plane')
	head['PLANEY_E'] = (error_y_new, 'scale factor error, error on the scale factor')
	head['PLN_CR1'] = (X_mid_h, 'reference centre pixel for plane in axis 1')
	head['PLN_CR2'] = (Y_mid_h, 'reference centre pixel for plane in axis 2')

	## Center pixel of the CSR
	

except ValueError:
	print ("ValueError: The two image are not same size.")
	exit(0)


plt.show()
print ()
print ('Done Plotting.')

# saving the files and the values 
if savefits == True:
	print ('\nSaving fits file for corrected data.')
	os.system("rm {0}".format(output_file))
	fits.writeto('{0}'.format(output_file), herschel_native, head)

if saveimage == True:	
	print ('\nSaving png file.')
	plt.savefig('{0}'.format(output_combine))

if savecc == True:
	print ('\nSaving colour correction fits.')
	os.system("rm {0}".format(cc_map))
	fits.writeto(cc_map, colour_corrections, head)	


plt.figure(figsize = (9,4))
plt.subplot(121)
residual_z = (planck-(herschel+b+z_plane)/a)
residuals_z = residual_z[np.isnan(residual_z)==False]

totalpixel = np.count_nonzero(~np.isnan(herschel)) 
bin_num = int(totalpixel/10)
n1, bins1, patches = plt.hist(residuals_z.ravel(), bins=bin_num, color='tomato',density=True, alpha = 0.7,align='right') 
popt1, pcov1 = curve_fit(f, bins1[1:], n1)
y1 = f(bins1[1:], *popt1)
plt.plot(bins1[1:], y1, c='k', linewidth=1, label = '{2}: $\mu = ${0}, $\sigma = ${1}'.format(round(popt1[1],4), round(popt1[2],4), wavelength))
plt.legend()
plt.title('With Plane + Residual cutoff')

plt.subplot(122)
plt.imshow(residual_z, origin = 'lower', vmin = maxp*res_lim*-1, vmax = maxp*res_lim, cmap = turbo_cmap); plt.colorbar()
plt.title('{0} at Wavelength {1}'.format(name, wavelength) + r'$\rm{\mu m}$')


if saveimage == True:	
	print ('\nSaving png file.')
	plt.savefig('{0}'.format(output_combine[:-4]+'.after.png'))


# print all the results
print ('FINAL RESULT')
print ('scale',a)
print ('zero-point:', b)
print ('error scale:',error_x)
print ('error zero-point:', error_y)


try:
	herschel_offset = head['H_OFFSET']
	print ('herschel derived offset'),head['H_OFFSET'] 
except:
	herschel_offset = float('Nan')
	print ('There is no offset provided by Herschel')

try:
	herschel_offset_err = head['H_OFF_E']
	print ('herschel derived offset error'),head['H_OFF_E']
except:
	herschel_offset_err = float('Nan')
	print ('There is no offset error provided by Herschel')


result = [str(planck_name), str(name), str(wavelength), 
		  str(a), str(error_a), str(b), str(error_b), str(herschel_offset), 
		  str(herschel_offset_err), str(scale_x), str(error_x_new), 
		  str(scale_y), str(error_y_new), str(X_mid_h), str(Y_mid_h) ]

result.append(str(np.nanmean(x_meas1)))
result.append(str(np.median(x_meas1)))
result.append(str(np.median(y_meas1)))
result.append(str(mean_before_crop))	
result.append(str(median_before_crop))
result.append(str(median_before_crop_her))


Result = {'Planck': [planck_name],
 		  'Herschel': [name], 
 		  'Wavelength': [wavelength],
 		  'Scale': [a],
 		  'Scale Error': [error_a],
 		  'Zero-point': [b],
 		  'Zero-point Error': [error_b],
 		  'SPIRE correction': [herschel_offset],
 		  'SPIRE Error': [herschel_offset_err],
 		  'Plane Scale X': [scale_x],
 		  'Scale X Error': [error_x_new],
 		  'Plane Scale Y': [scale_y],
 		  'Scale Y Error': [error_y_new],
 		  'Plane centre X': [X_mid_h],
 		  'Plane centre Y': [Y_mid_h],
 		  'Planck Median': [np.nanmedian(x_meas1)],
 		  'Herschel Median': [np.nanmedian(y_meas1)],
 		  'Planck Median Before': [median_before_crop],
 		  'Herschel Median Before': [median_before_crop_her]}

print ()
print (Result)

Column = ['Planck','Herschel', 'Wavelength','Scale','Scale Error',
		  'Zero-point','Zero-point Error', 'SPIRE correction','SPIRE Error',
		  'Plane Scale X', 'Scale X Error', 'Plane Scale Y', 'Scale Y Error', 
		  'Plane centre X', 'Plane centre Y',
		  'Planck Median' ,'Herschel Median' ,'Planck Median Before' ,'Herschel Median Before']

if write_offset == True:
	if os.path.exists(offset_csv):
	    print ("File exist")
	    header = False
	    #with open(offset_csv, 'a') as f:
		#	f.writelines(','.join(Column)+'\n')
		#	f.close()
	else:
	    print ("File not exist")
	    header = True	

	#print ("\nwriting to ")
	#print (offset_csv)
	#with open(offset_csv, 'a') as f:
	#	f.writelines('\t'.join(result)+'\n')
	#	f.close()

	df = pd.DataFrame(Result)
	mode = 'w' if header else 'a'
	df.to_csv(offset_csv, encoding='utf-8', mode=mode, header=header, index=False, columns=Column)
		

print ()
# End timer
end_timer = time.time()
# Print time to run
print ('Finished in ', str(datetime.timedelta(seconds=end_timer-start_timer)))
print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print ()

