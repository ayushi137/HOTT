from __future__ import print_function

from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rcParams

from turbo_colormap import turbo_colormap_data
from turbo_colormap import turbo_cmap
from matplotlib.colors import ListedColormap
from scipy.optimize import leastsq
from scipy.stats import chi2
import required_functions as rf

import sys
import os

import pandas as pd 


import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

# alternative to tight layout
rcParams.update({'figure.autolayout': True})
plt.ion()

import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# Equation for Gaussian
def f(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))

#fit function
def peval(x, p):
    return p[0] * np.exp(-(x - p[1])**2.0 / (2 * p[2]**2))

def residuals (p,y,x, peval):
    return (y) - peval(x,p)    


def getoffset(fieldname):
	df = pd.read_csv(offset_file,delimiter=',')
	y = df.index[df['Field']==fieldname]
	print ('The index for ', fieldname, ' is ', y[0])
	offset_values = np.array([df.iloc[y[0]]['Offset 160'],df.iloc[y[0]]['Offset 250'],df.iloc[y[0]]['Offset 350'],df.iloc[y[0]]['Offset 500']])
	return offset_values


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

	hist, bin_val = np.histogram(no_nan, bin_num)
	
	i = 0
	j = 0
	
	all_hist = hist
	all_bin = bin_val
	while i < np.sum(hist)*(low_percentage/100.):
		i += all_hist[j]
		j += 1
		#print (i,j)


	k = j 
	#print (j, k)

	return bin_val, k 



field = sys.argv[1]

save_files = True
save_files_offset = True
mask_map_used = True
high_crop = False
low_crop = False
percentage = 15
low_percentage = 50
tickfontsize = 15
fontsize = 15

if save_files == True:
	date = sys.argv[2]

#result_folder = '/Users/ayushisingh/Documents/Research/colden_project/'
#folder = '/Volumes/Storm/Research_has_all_data/colden_project/'
folder = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project/'
result_folder = folder


#gain_value =np.array([1.008, 0.994, 1.0, 1.0])
gain_value =np.array([1.0, 1.0, 1.0, 1.0])

initial_run = False
suffix = '_1_v2'
suffix2 = '_0_v2'
date_old = '2022-10-27'
offset_file = result_folder+'/offset_values/Offset_and_mu_{1}{0}.csv'.format(suffix2, date_old)

print ('Is it an initial run:', initial_run)
print ('suffix:',suffix)
if initial_run == False:
	print ('offset file:',offset_file)

print ()
grid = fits.getdata(folder+'grid/{0}_500.fits'.format(field))
grid[grid<0.9] = float('nan')
grid[grid>=0.9] = 1.


#### Load maps
w = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_model_160_orig.fits'.format(field,suffix))
x = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_model_250_orig.fits'.format(field,suffix))
y = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_model_350_orig.fits'.format(field,suffix))
z = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_model_500_orig.fits'.format(field,suffix))
chisq_map = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_chisquare_orig.fits'.format(field,suffix))

werr = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_error_160_orig.fits'.format(field,suffix))
xerr = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_error_250_orig.fits'.format(field,suffix))
yerr = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_error_350_orig.fits'.format(field,suffix))
zerr = fits.getdata(folder+'tau_temp/{0}{1}/{0}_hott_error_500_orig.fits'.format(field,suffix))

w = w/grid
x = x/grid
y = y/grid
z = z/grid
chisq_map = chisq_map/grid

W = fits.getdata(folder+'herschel_data/{0}/{0}-160.offset.conv.resamp.fits'.format(field))
X = fits.getdata(folder+'herschel_data/{0}/{0}-250.offset.conv.resamp.fits'.format(field))
Y = fits.getdata(folder+'herschel_data/{0}/{0}-350.offset.conv.resamp.fits'.format(field))
Z = fits.getdata(folder+'herschel_data/{0}/{0}-500.offset.fits'.format(field))

if mask_map_used == True:	
	try:
		maskW = fits.getdata(folder+'offset_mask/{0}-160.masked_mask.resamp.fits'.format(field))
		maskX = fits.getdata(folder+'offset_mask/{0}-250.masked_mask.resamp.fits'.format(field))
		maskY = fits.getdata(folder+'offset_mask/{0}-350.masked_mask.resamp.fits'.format(field))
		maskZ = fits.getdata(folder+'offset_mask/{0}-500.masked_mask.resamp.fits'.format(field))

	except:
		rf.regrid(folder+'offset_mask/{0}-160.masked_mask.fits'.format(field), folder+'herschel_data/{0}/{0}-160.offset.conv.resamp.fits'.format(field))
		rf.regrid(folder+'offset_mask/{0}-250.masked_mask.fits'.format(field), folder+'herschel_data/{0}/{0}-250.offset.conv.resamp.fits'.format(field))
		rf.regrid(folder+'offset_mask/{0}-350.masked_mask.fits'.format(field), folder+'herschel_data/{0}/{0}-350.offset.conv.resamp.fits'.format(field))
		rf.regrid(folder+'offset_mask/{0}-500.masked_mask.fits'.format(field), folder+'herschel_data/{0}/{0}-500.offset.fits'.format(field))

		maskW = fits.getdata(folder+'offset_mask/{0}-160.masked_mask.resamp.fits'.format(field))
		maskX = fits.getdata(folder+'offset_mask/{0}-250.masked_mask.resamp.fits'.format(field))
		maskY = fits.getdata(folder+'offset_mask/{0}-350.masked_mask.resamp.fits'.format(field))
		maskZ = fits.getdata(folder+'offset_mask/{0}-500.masked_mask.resamp.fits'.format(field))


print (folder+'herschel_data/{0}/{0}-160.offset.conv.resamp.fits'.format(field))

if initial_run == False:
	offset_values = getoffset(field)
	print ('offset_values:', offset_values)
	W = W*gain_value[0] - offset_values[0] 
	X = X*gain_value[1] - offset_values[1]
	Y = Y*gain_value[2] - offset_values[2]
	Z = Z*gain_value[3] - offset_values[3] 
else:
	print ('This is the initial run.')	

chi_mask = chisq_map[np.isnan(X)] = float('nan')
Whead = fits.getheader(folder+'herschel_data/{0}/{0}-160.offset.fits'.format(field))
Xhead = fits.getheader(folder+'herschel_data/{0}/{0}-250.offset.fits'.format(field))
Yhead = fits.getheader(folder+'herschel_data/{0}/{0}-350.offset.fits'.format(field))
Zhead = fits.getheader(folder+'herschel_data/{0}/{0}-500.offset.fits'.format(field))

Offset_W = Whead['OFFSET']
Offset_X = Xhead['OFFSET']
Offset_Y = Yhead['OFFSET']
Offset_Z = Zhead['OFFSET']

HerX = Xhead['H_OFFSET']
HerY = Yhead['H_OFFSET']
HerZ = Zhead['H_OFFSET']

Offset_X = Offset_X+HerX
Offset_Y = Offset_Y+HerY
Offset_Z = Offset_Z+HerZ


mask_chisq = np.isnan(chisq_map) == False
chisq = chisq_map[mask_chisq]
chi_mean = np.nanmean(chisq)
chi_median = np.nanmedian(chisq)

cmap = turbo_cmap
colors = cmap(np.linspace(0, 1, (4)))
c = colors[0:4]
#c = c[::-1]

####################################################################

LW = 2
legend_size = 8

print ()
print ('Getting Data-model')

plt.figure(figsize = (14,8))
plt.subplot(231)
plt.title(field, fontsize = fontsize)
rw = (W-w)
rx = (X-x)
ry = (Y-y)
rz = (Z-z)

if high_crop == True:

	cuttoff_map = np.copy(w)
	bin_val, k = high_cropping_fun(cuttoff_map, percentage)
	rw[cuttoff_map > bin_val[k]] = float('nan')

	cuttoff_map = np.copy(x)
	bin_val, k = high_cropping_fun(cuttoff_map, percentage)
	rx[cuttoff_map > bin_val[k]] = float('nan')

	cuttoff_map = np.copy(y)
	bin_val, k = high_cropping_fun(cuttoff_map, percentage)
	ry[cuttoff_map > bin_val[k]] = float('nan')

	cuttoff_map = np.copy(z)
	bin_val, k = high_cropping_fun(cuttoff_map, percentage)
	rz[cuttoff_map > bin_val[k]] = float('nan')

	print ('Finish removing the higher points.')


if mask_map_used == True:	
	print ('mask_map_used is True and high_crop is False')

	maskW[maskW<0.5] = float('nan')
	maskW[maskW>=0.5] = 1.0
	maskX[maskX<0.5] = float('nan')
	maskX[maskX>=0.5] = 1.0
	maskY[maskY<0.5] = float('nan')
	maskY[maskY>=0.5] = 1.0
	maskZ[maskZ<0.5] = float('nan')
	maskZ[maskZ>=0.5] = 1.0

	rw =rw/maskW
	rx =rx/maskX
	ry =ry/maskY
	rz =rz/maskZ


mask_pixel= (np.isnan(rw)==False)& (np.isnan(rx)==False)& (np.isnan(ry)==False)& (np.isnan(rz)==False)

rW = (rw[mask_pixel])
rX = (rx[mask_pixel])
rY = (ry[mask_pixel])
rZ = (rz[mask_pixel])


totalpixel = np.count_nonzero(~np.isnan(rW))
bin_num = int(totalpixel/10)
print ('bins:', bin_num)
#plt.subplot(2,2,1)
n1, bins1, patches = plt.hist(rW, bins=bin_num, range= (-60.0, 60.0),color=c[0],density=True, alpha = 0.3)
n2, bins2, patches = plt.hist(rX, bins=bin_num, range= (-10.0, 10.0),color=c[1],density=True, alpha = 0.3)
n3, bins3, patches = plt.hist(rY, bins=bin_num, range= (-10.0, 10.0),color=c[2],density=True, alpha = 0.3)
n4, bins4, patches = plt.hist(rZ, bins=bin_num, range= (-10.0, 10.0),color=c[3],density=True, alpha = 0.3)

plt.plot(0,0,'w.', label = r'  $\lambda$:        $\mu$,        $\sigma$')
p0 = [0.0, np.median(rW), np.std(rW)]
popt1 = leastsq(residuals,p0,args=(n1,bins1[:-1], peval))
y1 = peval(bins1[:-1],popt1[0])
plt.plot(bins1[:-1], y1, c=c[0], linewidth=LW, label = '160: {0},  {1}'.format(round(popt1[0][1],4), round(popt1[0][2],2)))

p0 = [0.0, np.median(rX), np.std(rX)]
popt2 = leastsq(residuals,p0,args=(n2,bins2[:-1], peval))
y2 = peval(bins2[:-1],popt2[0])
plt.plot(bins2[:-1], y2, c=c[1], linewidth=LW, label = '250: {0},  {1}'.format(round(popt2[0][1],4), round(popt2[0][2],2)))

p0 = [0.0, np.median(rY), np.std(rY)]
popt3 = leastsq(residuals,p0,args=(n3,bins3[:-1], peval))
y3 = peval(bins3[:-1],popt3[0])
plt.plot(bins3[:-1], y3, c=c[2], linewidth=LW, label = '350: {0},  {1}'.format(round(popt3[0][1],4), round(popt3[0][2],2)))

p0 = [0.0, np.median(rZ), np.std(rZ)]
popt4 = leastsq(residuals,p0,args=(n4,bins4[:-1], peval))
y4 = peval(bins4[:-1],popt4[0])
plt.plot(bins4[:-1], y4, c=c[3], linewidth=LW, label = '500: {0},  {1}'.format(round(popt4[0][1],4), round(popt4[0][2],2)))

#plt.title('{0}'.format(field))
plt.legend(loc=2,ncol=1,prop={'size':legend_size})
#plt.xlim(-5e20, 5e20)
plt.xlim(-3,3)
plt.xlabel(r'$I_{\nu} - I_{\nu,m}$', fontsize = fontsize)
#plt.xlabel(r'[MJy/sr]', fontsize = fontsize)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)

DMD =  [str(popt1[0][1]), str(popt2[0][1]), str(popt3[0][1]), str(popt4[0][1])]
DMD_sigma =  [str(popt1[0][2]), str(popt2[0][2]), str(popt3[0][2]), str(popt4[0][2])]
print ('Offsets are:',DMD)





####################################################################
print ()
print ('Getting Data-model/Data')
plt.subplot(232)

rw = (W-w)/W
rx = (X-x)/X
ry = (Y-y)/Y
rz = (Z-z)/Z


if low_crop == True:
	
	low_cuttoff_map = np.copy(w)
	bin_val, k = low_cropping_fun(low_cuttoff_map, low_percentage)
	rw[low_cuttoff_map < bin_val[k]] = float('nan')

	low_cuttoff_map = np.copy(x)
	bin_val, k = low_cropping_fun(low_cuttoff_map, low_percentage)
	rx[low_cuttoff_map < bin_val[k]] = float('nan')

	low_cuttoff_map = np.copy(y)
	bin_val, k = low_cropping_fun(low_cuttoff_map, low_percentage)
	ry[low_cuttoff_map < bin_val[k]] = float('nan')

	low_cuttoff_map = np.copy(z)
	bin_val, k = low_cropping_fun(low_cuttoff_map, low_percentage)
	rz[low_cuttoff_map < bin_val[k]] = float('nan')

	print ('Finish removing the lower points.')



mask_pixel= (np.isnan(rw)==False)& (np.isnan(rx)==False)& (np.isnan(ry)==False)& (np.isnan(rz)==False)

rW = (rw[mask_pixel])
rX = (rx[mask_pixel])
rY = (ry[mask_pixel])
rZ = (rz[mask_pixel])

totalpixel = np.count_nonzero(~np.isnan(rW))
bin_num = int(totalpixel/10)#200#int(totalpixel/2)
print ('bins:', bin_num)
#plt.subplot(2,2,1)
n1, bins1, patches = plt.hist(rW, bins=bin_num, range= (np.min(rW), np.max(rW)),color=c[0],density=True, alpha = 0.3)
n2, bins2, patches = plt.hist(rX, bins=bin_num, range= (np.min(rX), np.max(rX)),color=c[1],density=True, alpha = 0.3)
n3, bins3, patches = plt.hist(rY, bins=bin_num, range= (np.min(rY), np.max(rY)),color=c[2],density=True, alpha = 0.3)
n4, bins4, patches = plt.hist(rZ, bins=bin_num, range= (np.min(rZ), np.max(rZ)),color=c[3],density=True, alpha = 0.3)

plt.plot(0,0,'w.', label = r'  $\lambda$:        $\mu$,        $\sigma$')
plt.vlines(0, 0, np.max([n1, n2, n3, n4]), color= 'black',linestyle = '--') 
p0 = [0.0, np.median(rW), np.std(rW)]
popt1 = leastsq(residuals,p0,args=(n1,bins1[:-1], peval))
y1 = peval(bins1[:-1],popt1[0])
plt.plot(bins1[:-1], y1, c=c[0], linewidth=LW, label = '160: {0},  {1}'.format(round(popt1[0][1],3), round(popt1[0][2],3)))

p0 = [0.0, np.median(rX), np.std(rX)]
popt2 = leastsq(residuals,p0,args=(n2,bins2[:-1], peval))
y2 = peval(bins2[:-1],popt2[0])
plt.plot(bins2[:-1], y2, c=c[1], linewidth=LW, label = '250: {0},  {1}'.format(round(popt2[0][1],3), round(popt2[0][2],3)))

p0 = [0.0, np.median(rY), np.std(rY)]
popt3 = leastsq(residuals,p0,args=(n3,bins3[:-1], peval))
y3 = peval(bins3[:-1],popt3[0])
plt.plot(bins3[:-1], y3, c=c[2], linewidth=LW, label = '350: {0},  {1}'.format(round(popt3[0][1],3), round(popt3[0][2],3)))

p0 = [0.0, np.median(rZ), np.std(rZ)]
popt4 = leastsq(residuals,p0,args=(n4,bins4[:-1], peval))
y4 = peval(bins4[:-1],popt4[0])
plt.plot(bins4[:-1], y4, c=c[3], linewidth=LW, label = '500: {0},  {1}'.format(round(popt4[0][1],3), round(popt4[0][2],3)))

#plt.title('{0}'.format(field))
plt.legend(loc=2,ncol=1,prop={'size':legend_size})
#plt.xlim(-5e20, 5e20)
plt.xlim(-0.05, 0.05)
plt.xlabel(r'$(I_{\nu} - I_{\nu,m})/I_{\nu}$', fontsize = fontsize)
#plt.xlabel(r'[MJy/sr]', fontsize = fontsize)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)

D_M_over_D = [popt1[0][1], popt2[0][1], popt3[0][1], popt4[0][1]]
print ('Offsets are:',D_M_over_D)


####################################################################


print ()
print ('Getting Data-model/error')

plt.subplot(233)
rw = (W-w)/werr
rx = (X-x)/xerr
ry = (Y-y)/yerr
rz = (Z-z)/zerr

mask_pixel= (np.isnan(rw)==False)& (np.isnan(rx)==False)& (np.isnan(ry)==False)& (np.isnan(rz)==False)

rW = (rw[mask_pixel])
rX = (rx[mask_pixel])
rY = (ry[mask_pixel])
rZ = (rz[mask_pixel])

totalpixel = np.count_nonzero(~np.isnan(rW))
bin_num = int(totalpixel/100)

print ('bins:', bin_num)
#plt.subplot(2,2,1)
n1, bins1, patches = plt.hist(rW, bins=bin_num, range= (-5.0, 5.0),color=c[0],density=True, alpha = 0.3)
n2, bins2, patches = plt.hist(rX, bins=bin_num, range= (-5.0, 5.0),color=c[1],density=True, alpha = 0.3)
n3, bins3, patches = plt.hist(rY, bins=bin_num, range= (-5.0, 5.0),color=c[2],density=True, alpha = 0.3)
n4, bins4, patches = plt.hist(rZ, bins=bin_num, range= (-5.0, 5.0),color=c[3],density=True, alpha = 0.3)

plt.plot(0,0,'w.', label = r'  $\lambda$:        $\mu$,        $\sigma$')
p0 = [0.0, np.median(rW), np.std(rW)]
popt1 = leastsq(residuals,p0,args=(n1,bins1[:-1], peval))
y1 = peval(bins1[:-1],popt1[0])
plt.plot(bins1[:-1], y1, c=c[0], linewidth=LW, label = '160: {0},  {1}'.format(round(popt1[0][1],3), round(popt1[0][2],2)))

p0 = [0.0, np.median(rX), np.std(rX)]
popt2 = leastsq(residuals,p0,args=(n2,bins2[:-1], peval))
y2 = peval(bins2[:-1],popt2[0])
plt.plot(bins2[:-1], y2, c=c[1], linewidth=LW, label = '250: {0},  {1}'.format(round(popt2[0][1],3), round(popt2[0][2],2)))

p0 = [0.0, np.median(rY), np.std(rY)]
popt3 = leastsq(residuals,p0,args=(n3,bins3[:-1], peval))
y3 = peval(bins3[:-1],popt3[0])
plt.plot(bins3[:-1], y3, c=c[2], linewidth=LW, label = '350: {0},  {1}'.format(round(popt3[0][1],3), round(popt3[0][2],2)))

p0 = [0.0, np.median(rZ), np.std(rZ)]
popt4 = leastsq(residuals,p0,args=(n4,bins4[:-1], peval))
y4 = peval(bins4[:-1],popt4[0])
plt.plot(bins4[:-1], y4, c=c[3], linewidth=LW, label = '500: {0},  {1}'.format(round(popt4[0][1],3), round(popt4[0][2],2)))

#plt.title('{0}'.format(field))
plt.legend(loc=2,ncol=1,prop={'size':legend_size})
#plt.xlim(-5e20, 5e20)
plt.xlim(-3.0, 3.0)

#plt.title(r'$(I_{\nu} - I_{\nu,m})/\sigma_{\nu}$')
plt.xlabel(r'$(I_{\nu} - I_{\nu,m})/\sigma_{\nu}$', fontsize = fontsize)
#plt.xlabel(r'[MJy/sr]', fontsize = fontsize)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)

D_M_over_D_sigma = [popt1[0][2], popt2[0][2], popt3[0][2], popt4[0][2]]
print ('Offsets are:',D_M_over_D_sigma)

####################################################################
print ('Getting second row')

cuttoff_map = np.copy(chisq_map)
mastermask = np.sum(np.array([maskW, maskX, maskY, maskZ]), axis = 0)

mastermask_mod = np.copy(mastermask)/grid
mastermask_mod[mastermask_mod<=3.5] = float('nan')
mastermask_mod[mastermask_mod>3.5] = 1.0
hist_map_mod = chisq_map/mastermask_mod
mask_chisq_new = (np.isnan(hist_map_mod) == False)
chisq_masked=hist_map_mod[mask_chisq_new]


mastermask_mod2 = np.copy(mastermask)
mastermask_mod2[np.isnan(mastermask_mod2)==True] = 1.0
mastermask_mod2[mastermask_mod2>=3.5] = float('nan')
mastermask_mod2 = mastermask_mod2/grid
hist_map_mod2 = chisq_map/mastermask_mod2
mask_chisq_new2 = (np.isnan(hist_map_mod2) == False)
chisq_masked2=hist_map_mod2[mask_chisq_new2]

###################################################################
print ('Getting Chi_map')

bin_val, k = high_cropping_fun(cuttoff_map, 10)
cuttoff_map[cuttoff_map > bin_val[k]] = float('nan')

plt.subplot(234)
im = plt.imshow(chisq_map, origin = 'lower', cmap=cmap, vmax = np.nanmax(cuttoff_map)*0.6); 
cbar = plt.colorbar(im, pad=0.01)
cbar.ax.tick_params(labelsize = tickfontsize-2)
plt.title(r'$\chi^2$'.format(field), fontsize = fontsize)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)
plt.tick_params(labelsize = tickfontsize)

####################################################################

plt.subplot(235)
N = 256
vals2 = np.ones((N, 4))
vals2[:, 0] = np.linspace(c[2][0], 1, N)
vals2[:, 1] = np.linspace(c[2][1], 1, N)
vals2[:, 2] = np.linspace(c[2][2], 1, N)
newcmp2 = ListedColormap(vals2)

N = 256
vals1 = np.ones((N, 4))
vals1[:, 0] = np.linspace(c[1][0], 1, N)
vals1[:, 1] = np.linspace(c[1][1], 1, N)
vals1[:, 2] = np.linspace(c[1][2], 1, N)
newcmp1 = ListedColormap(vals1)

N = 256
vals3 = np.ones((N, 4))
vals3[:, 0] = np.linspace(c[3][0], 1, N)
vals3[:, 1] = np.linspace(c[3][1], 1, N)
vals3[:, 2] = np.linspace(c[3][2], 1, N)
newcmp3 = ListedColormap(vals3)

N = 256
vals4 = np.ones((N, 4))
vals4[:, 0] = np.linspace(c[1][0], c[2][0], N)
vals4[:, 1] = np.linspace(c[1][1], c[2][1], N)
vals4[:, 2] = np.linspace(c[1][2], c[2][2], N)
newcmp4 = ListedColormap(vals4)

flat_mask = np.copy(mastermask)*0+1
flat_mask[np.isnan(flat_mask)== True] = 0
flat_mask2 = np.copy(flat_mask)
flat_mask[flat_mask == 1] = float('nan')
flat_mask = flat_mask/grid

#plt.imshow(np.log10(X/grid), origin='lower', cmap = "Greys", vmin = 0, vmax = np.log10(np.nanmax(X/grid)*0.6)); plt.colorbar()  
#plt.imshow((grid), origin='lower', cmap = newcmp3, alpha = 0.3);# plt.colorbar() 
#plt.imshow(flat_mask2/grid, origin='lower', cmap = newcmp4, alpha = 0.3)                                                                                 

im = plt.imshow((X/grid), origin='lower', cmap = "Greys", vmin = 0, vmax = np.nanmax(X/grid)*0.3); 
cbar = plt.colorbar(im, pad=0.01)
cbar.set_label(r'MJy sr$^{-1}$',size= fontsize,rotation=270,labelpad=20)
cbar.ax.tick_params(labelsize = tickfontsize-2)

plt.imshow((grid), origin='lower', cmap = newcmp3, alpha = 0.2);# plt.colorbar() 
plt.imshow(flat_mask2/grid, origin='lower', cmap = newcmp4, alpha = 0.3)         
plt.contour( flat_mask2/grid, [0,1], cmap = newcmp4, alpha = 0.4)


#plt.imshow((grid), origin='lower', cmap = newcmp3, alpha = 0.7);# plt.colorbar() 
#plt.imshow((grid), origin='lower', cmap = newcmp2, alpha = 0.2);# plt.colorbar() 
#plt.imshow((flat_mask), origin='lower', cmap = newcmp1, alpha = 0.2);# plt.colorbar() 

plt.title('Union mask', fontsize = fontsize)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)
plt.tick_params(labelsize = tickfontsize)

####################################################################
print ()
print ('Getting Chi_sq')
#plt.figure(figsize = (7,5))
totalpixel = np.count_nonzero(~np.isnan(chisq))
bin_num = int(totalpixel/700)
print ('bins:', bin_num)
plt.subplot(236)



# 1,0,3

range_max = 5
range_min = -3
deltax = range_max/bin_num


n0, bins0 = np.histogram(np.log10(chisq), bins=bin_num, range=(range_min,range_max))


n0_0, bins0_0, patches = plt.hist(np.log10(chisq), bins=bin_num, range=(range_min,range_max),color=c[3], alpha = 0.7, density = True, label = r'Total')


norm_value = (n0/n0_0)
norm_value = norm_value[0]

counts1, binsnp1 = np.histogram(np.log10(chisq_masked), bins=bin_num,range=(range_min,range_max))

n1, bins1, patches1 = plt.hist(binsnp1[:-1], binsnp1, weights=counts1/norm_value , color=c[2], alpha = 0.7, label = r'Low')


counts2, binsnp2 = np.histogram(np.log10(chisq_masked2), bins=bin_num,range=(range_min,range_max))

n2, bins2, patches2 = plt.hist(binsnp2[:-1], binsnp2, weights=counts2/norm_value,color=c[1], alpha = 0.7, label = r'High')





p0 = [0.0, 1.2, 0.5]

popt0 = leastsq(residuals,p0,args=(n0,bins0[:-1], peval))
#y0 = peval(bins0[:-1],popt0[0])
#plt.plot(bins0[:-1], y0, c=c[2], linewidth=1)

#x_vals= np.arange(0, np.nanmax(chisq), 0.01)
x_vals= np.arange(range_min, range_max, 0.1)
x_vals = 10**x_vals
chi2_0 = chi2.stats(2, moments='mvsk')[0]
chi2_pdf = chi2.pdf(x_vals, 2)
plt.legend(prop={'size':legend_size+2})
plt.plot(np.log10(x_vals), np.log(10)*x_vals*chi2_pdf, c='k', linewidth = LW, label = r'$\chi^2$ at k = 2')
#plt.xlabel(r'$\chi^2$'.format(field), fontsize = fontsize)
plt.xlabel(r'$\rm{log}(\chi^2)$', fontsize = fontsize)
plt.xlim(-3, 2)
plt.xticks( fontsize=tickfontsize)
plt.yticks( fontsize=tickfontsize)


####################################################################


plt.tight_layout()


result = {'Field': [field],
		  'Zero-point 160': [Offset_W],
		  'Zero-point 250': [Offset_X],
		  'Zero-point 350': [Offset_Y],
		  'Zero-point 500': [Offset_Z],
		  'Offset 160': [DMD[0]],
		  'Offset 250': [DMD[1]],
		  'Offset 350': [DMD[2]],
		  'Offset 500': [DMD[3]],
		  'Offset 160 sigma': [DMD_sigma[0]],
		  'Offset 250 sigma': [DMD_sigma[1]],
		  'Offset 350 sigma': [DMD_sigma[2]],
		  'Offset 500 sigma': [DMD_sigma[3]],
		  'Gain 160': [D_M_over_D[0]],
		  'Gain 250': [D_M_over_D[1]],
		  'Gain 350': [D_M_over_D[2]],
		  'Gain 500': [D_M_over_D[3]],
		  'Gain 160 sigma': [D_M_over_D_sigma[0]],
		  'Gain 250 sigma': [D_M_over_D_sigma[1]],
		  'Gain 350 sigma': [D_M_over_D_sigma[2]],
		  'Gain 500 sigma': [D_M_over_D_sigma[3]]}

Column = ['Field','Zero-point 160','Zero-point 250','Zero-point 350','Zero-point 500',
'Offset 160','Offset 250','Offset 350','Offset 500', 'Offset 160 sigma','Offset 250 sigma','Offset 350 sigma','Offset 500 sigma',
'Gain 160','Gain 250','Gain 350','Gain 500', 'Gain 160 sigma','Gain 250 sigma','Gain 350 sigma','Gain 500 sigma']
	

if save_files == True:
	os.system('mkdir {1}/histograms/histograms_{0}'.format(date,result_folder))
	
	plt.savefig(result_folder+'histograms/histograms_{1}/{0}{2}.png'.format(field,date,suffix))
	plt.savefig(result_folder+'histograms/histograms_{1}/{0}{2}.pdf'.format(field,date,suffix))

	offset_csv = result_folder+'offset_values/Offset_and_mu_'+date+suffix+'.csv'

	if os.path.exists(offset_csv):
	    print ("File exist")
	    header = False
	   	
	else:
	    print ("File not exist")
	    header = True


	df = pd.DataFrame(result)
	mode = 'w' if header else 'a'
	df.to_csv(offset_csv, encoding='utf-8', mode=mode, header=header, index=False, columns=Column)
	


