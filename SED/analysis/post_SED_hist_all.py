from __future__ import print_function

from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from sys import exit
from turbo_colormap import turbo_colormap_data
from matplotlib.colors import ListedColormap

import sys
import os
import pandas as pd 
import time, datetime
    
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
# alternative to tight layout
rcParams.update({'figure.autolayout': True})
plt.ion()

# Equation for Gaussian
def f(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))

individual = False
histOfhist = False
combined = False
full = False
paper = True

Column = ['Field','Zero-point 160','Zero-point 250','Zero-point 350','Zero-point 500',
'Offset 160 sigma','Offset 250 sigma','Offset 350 sigma','Offset 500 sigma', 'Offset 160','Offset 250','Offset 350','Offset 500',
'Gain 160','Gain 250','Gain 350','Gain 500', 'Gain 160 sigma','Gain 250 sigma','Gain 350 sigma','Gain 500 sigma']

date = str((datetime.date.today()))
#folder = '/Users/ayushisingh/Documents/Research/colden_project/offset_values/'
folder = '/Volumes/Storm/Research_has_all_data/colden_project/offset_values/'
#result = '/Volumes/Storm/Research_has_all_data/colden_project/histograms/results_combined_hist_'+date+'_1_9/'
result = '/Volumes/Storm/Research_has_all_data/colden_project/histograms/results_combined_hist_'+date+'_1_1/'
print (result)

os.system('mkdir '+result)

suffix = ''#'_mask2'

#offset_file = folder+'Offset_and_mu_2021-06-09_1_9.csv'#09-11_1_8.csv'
#offset_file = folder+'Offset_and_mu_2022-03-21_1_1.csv'#09-11_1_8.csv'
offset_file = folder+'Offset_and_mu_2022-06-03_1_1.csv'#09-11_1_8.csv'


df = pd.read_csv(offset_file,delimiter=',')

name = np.array(df.Field)
mu_sub_160 = np.array(df['Offset 160'])
mu_sub_250 = np.array(df['Offset 250'])
mu_sub_350 = np.array(df['Offset 350'])
mu_sub_500 = np.array(df['Offset 500'])
mu_sub_160_sigma = np.array(df['Offset 160 sigma'])
mu_sub_250_sigma = np.array(df['Offset 250 sigma'])
mu_sub_350_sigma = np.array(df['Offset 350 sigma'])
mu_sub_500_sigma = np.array(df['Offset 500 sigma'])
offset_160 = np.array(df['Zero-point 160'])
offset_250 = np.array(df['Zero-point 250'])
offset_350 = np.array(df['Zero-point 350'])
offset_500 = np.array(df['Zero-point 500'])
mu_160 = np.array(df['Gain 160'])
mu_250 = np.array(df['Gain 250'])
mu_350 = np.array(df['Gain 350'])
mu_500 = np.array(df['Gain 500'])
mu_160_sigma = np.array(df['Gain 160 sigma'])
mu_250_sigma = np.array(df['Gain 250 sigma'])
mu_350_sigma = np.array(df['Gain 350 sigma'])
mu_500_sigma = np.array(df['Gain 500 sigma'])


for i in range(len(name)):
    print (i, name[i], '\t\t', round(offset_160[i],2),round(offset_250[i],2),round(offset_350[i],2),round(offset_500[i],2)) 

big = np.array([18,38,54,55,56,57])#14,38])
ia = np.indices(mu_sub_160.shape)
small = np.setxor1d(ia, big)

big = big.astype('int64')
small = small.astype('int64') 

#plt.xlim(-0.2, 0.2)
#plt.title('{0}: Data-model/Data'.format(field))
#plt.xlabel(r'[MJy/sr]')

cmap = ListedColormap(turbo_colormap_data)#cm.inferno
colors = cmap(np.linspace(0, 1, (4)))
c = colors[0:4]

print ()
print (big)
print ()
for i in range(len(name)): 
    if i in big: 
        print (name[i])

if paper == True:
    
    range_val = 0.015
    BINS = 40
    plt.figure(figsize = (7,5))
    #plt.subplot2grid((2, 4), (0, 0))#plt.subplot(2,4,1)
    #n1, bins1, patches = plt.hist(mu_160[small], bins=30, range= (-0.05, 0.05), color=c[0],density=True, alpha = 0.7,align='right')
    #plt.subplot2grid((2, 4), (0, 1))
    #n2, bins2, patches = plt.hist(mu_250[small], bins=30, range= (-0.05, 0.05),color=c[1],density=True, alpha = 0.7,align='right')
    #plt.subplot2grid((2, 4), (1, 0))
    #n3, bins3, patches = plt.hist(mu_350[small], bins=30, range= (-0.05, 0.05),color=c[2],density=True, alpha = 0.7,align='right')
    #plt.subplot2grid((2, 4), (1, 1))
    #n4, bins4, patches = plt.hist(mu_500[small], bins=30, range= (-0.05, 0.05),color=c[3],density=True, alpha = 0.7,align='right')

    #plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=2)
    n1, bins1, patches = plt.hist(mu_160[small], bins=BINS, range= (-range_val, range_val), color=c[0], alpha = 0.7,align='right')
    #plt.subplot(2,2,2)
    n2, bins2, patches = plt.hist(mu_250[small], bins=BINS, range= (-range_val, range_val),color=c[1], alpha = 0.7,align='right')
    #plt.subplot(2,2,3)
    n3, bins3, patches = plt.hist(mu_350[small], bins=BINS, range= (-range_val, range_val),color=c[2], alpha = 0.5,align='right')
    #plt.subplot(2,2,4)
    n4, bins4, patches = plt.hist(mu_500[small], bins=BINS, range= (-range_val, range_val),color=c[3], alpha = 0.5,align='right')

    
    plt.plot(0,0,'w.', label = r'  $\lambda$:        $\mu$,          $\sigma$')
    #plt.plot(0,0,'w.', label = r'  $\lambda$:        mean,          std')
    popt1, pcov1 = curve_fit(f, bins1[1:], n1)
    y1 = f(bins1[1:], *popt1)
    plt.plot(bins1[1:], y1, c=c[0], linewidth=2, label = '160: {0},  {1}'.format(round(popt1[1],3), abs(round(popt1[2],3))))

    popt2, pcov2 = curve_fit(f, bins2[1:], n2)
    y2 = f(bins2[1:], *popt2)
    plt.plot(bins2[1:], y2, c=c[1], linewidth=2, label = '250: {0},  {1}'.format(round(popt2[1],3), abs(round(popt2[2],3))))

    popt3, pcov3 = curve_fit(f, bins3[1:], n3)
    y3 = f(bins3[1:], *popt3)
    plt.plot(bins3[1:], y3, c=c[2], linewidth=2, label = '350: {0},  {1}'.format(round(popt3[1],3), abs(round(popt3[2],3))))

    popt4, pcov4 = curve_fit(f, bins4[1:], n4)
    y4 = f(bins4[1:], *popt4)
    plt.plot(bins4[1:], y4, c=c[3], linewidth=2, label = '500: {0},  {1}'.format(round(popt4[1],3), abs(round(popt4[2],3))))
    
    plt.xlabel(r'$(I_{\nu} - I_{\nu,m})/I_{\nu}$',fontsize= 15)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize= 15)
    plt.legend(loc=2,ncol=1,prop={'size':12})
    #plt.xlim(-5e20, 5e20)
    plt.xlim(-range_val, range_val)
    plt.tick_params(labelsize = 15)
    plt.savefig(result+'All_Gain{0}_main.png'.format(suffix))
    #plt.savefig(main_folder+'histograms/{0}.fits'.format(field))

    print(np.nanmean(mu_160[small]))
    print(np.nanmean(mu_250[small]))
    print(np.nanmean(mu_350[small]))
    print(np.nanmean(mu_500[small]))
    print()
    print(np.nanmedian(mu_160[small]))
    print(np.nanmedian(mu_250[small]))
    print(np.nanmedian(mu_350[small]))
    print(np.nanmedian(mu_500[small]))

# to get the histograms of all the mu together or speperate. 
if histOfhist == True:
    plt.figure(figsize = (10,5))
    plt.subplot2grid((2, 4), (0, 0))#plt.subplot(2,4,1)
    n1, bins1, patches = plt.hist(mu_sub_160[small], bins=30, range= (-3.0, 3.0), color=c[0],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (0, 1))
    n2, bins2, patches = plt.hist(mu_sub_250[small], bins=30, range= (-3.0, 3.0),color=c[1],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (1, 0))
    n3, bins3, patches = plt.hist(mu_sub_350[small], bins=50, range= (-1.0, 1.0),color=c[2],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (1, 1))
    n4, bins4, patches = plt.hist(mu_sub_500[small], bins=50, range= (-1.0, 1.0),color=c[3],density=True, alpha = 0.7,align='right')

    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=2)
    n1, bins1, patches = plt.hist(mu_sub_160[small], bins=30, range= (-3.0, 3.0), color=c[0],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,2)
    n2, bins2, patches = plt.hist(mu_sub_250[small], bins=30, range= (-3.0, 3.0),color=c[1],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,3)
    n3, bins3, patches = plt.hist(mu_sub_350[small], bins=50, range= (-1.0, 1.0),color=c[2],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,4)
    n4, bins4, patches = plt.hist(mu_sub_500[small], bins=50, range= (-1.0, 1.0),color=c[3],density=True, alpha = 0.7,align='right')

    print ('160:',np.mean(mu_sub_160[small]), np.std(mu_sub_160[small]))
    print ('250:',np.mean(mu_sub_250[small]), np.std(mu_sub_250[small]))
    print ('350:',np.mean(mu_sub_350[small]), np.std(mu_sub_350[small]))
    print ('500:',np.mean(mu_sub_500[small]), np.std(mu_sub_500[small]))
    '''
    popt1, pcov1 = curve_fit(f, bins1[1:], n1)
    y1 = f(bins1[1:], *popt1)
    plt.plot(bins1[1:], y1, c=c[0], linewidth=1, label = '160: $\mu = ${0}, $\sigma = ${1}'.format(round(popt1[1],4), round(popt1[2],4)))

    popt2, pcov2 = curve_fit(f, bins2[1:], n2)
    y2 = f(bins2[1:], *popt2)
    plt.plot(bins2[1:], y2, c=c[1], linewidth=1, label = '250: $\mu = ${0}, $\sigma = ${1}'.format(round(popt2[1],4), round(popt2[2],4)))

    popt3, pcov3 = curve_fit(f, bins3[1:], n3)
    y3 = f(bins3[1:], *popt3)
    plt.plot(bins3[1:], y3, c=c[2], linewidth=1, label = '350: $\mu = ${0}, $\sigma = ${1}'.format(round(popt3[1],4), round(popt3[2],4)))

    popt4, pcov4 = curve_fit(f, bins4[1:], n4)
    y4 = f(bins4[1:], *popt4)
    plt.plot(bins4[1:], y4, c=c[3], linewidth=1, label = '500: $\mu = ${0}, $\sigma = ${1}'.format(round(popt4[1],4), round(popt4[2],4)))
    '''
    plt.title('All of them Data-Model')
    plt.legend(loc=2,ncol=1,prop={'size':10})
    #plt.xlim(-5e20, 5e20)
    plt.xlim(-5.0, 5.0)
    plt.savefig(result+'All_data-model{0}.png'.format(suffix))

    plt.figure(figsize = (10,5))
    plt.subplot2grid((2, 4), (0, 0))#plt.subplot(2,4,1)
    n1, bins1, patches = plt.hist(mu_160[small], bins=30, range= (-0.05, 0.05), color=c[0],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (0, 1))
    n2, bins2, patches = plt.hist(mu_250[small], bins=30, range= (-0.05, 0.05),color=c[1],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (1, 0))
    n3, bins3, patches = plt.hist(mu_350[small], bins=30, range= (-0.05, 0.05),color=c[2],density=True, alpha = 0.7,align='right')
    plt.subplot2grid((2, 4), (1, 1))
    n4, bins4, patches = plt.hist(mu_500[small], bins=30, range= (-0.05, 0.05),color=c[3],density=True, alpha = 0.7,align='right')

    plt.subplot2grid((2, 4), (0, 2), colspan=2, rowspan=2)
    n1, bins1, patches = plt.hist(mu_160[small], bins=30, range= (-0.05, 0.05), color=c[0],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,2)
    n2, bins2, patches = plt.hist(mu_250[small], bins=30, range= (-0.05, 0.05),color=c[1],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,3)
    n3, bins3, patches = plt.hist(mu_350[small], bins=30, range= (-0.05, 0.05),color=c[2],density=True, alpha = 0.7,align='right')
    #plt.subplot(2,2,4)
    n4, bins4, patches = plt.hist(mu_500[small], bins=30, range= (-0.05, 0.05),color=c[3],density=True, alpha = 0.7,align='right')

    
    popt1, pcov1 = curve_fit(f, bins1[1:], n1)
    y1 = f(bins1[1:], *popt1)
    plt.plot(bins1[1:], y1, c=c[0], linewidth=1, label = '160: $\mu = ${0}, $\sigma = ${1}'.format(round(popt1[1],4), round(popt1[2],4)))

    popt2, pcov2 = curve_fit(f, bins2[1:], n2)
    y2 = f(bins2[1:], *popt2)
    plt.plot(bins2[1:], y2, c=c[1], linewidth=1, label = '250: $\mu = ${0}, $\sigma = ${1}'.format(round(popt2[1],4), round(popt2[2],4)))

    popt3, pcov3 = curve_fit(f, bins3[1:], n3)
    y3 = f(bins3[1:], *popt3)
    plt.plot(bins3[1:], y3, c=c[2], linewidth=1, label = '350: $\mu = ${0}, $\sigma = ${1}'.format(round(popt3[1],4), round(popt3[2],4)))

    popt4, pcov4 = curve_fit(f, bins4[1:], n4)
    y4 = f(bins4[1:], *popt4)
    plt.plot(bins4[1:], y4, c=c[3], linewidth=1, label = '500: $\mu = ${0}, $\sigma = ${1}'.format(round(popt4[1],4), round(popt4[2],4)))
    
    plt.title('All of them Data-Model/Data')
    plt.legend(loc=2,ncol=1,prop={'size':10})
    #plt.xlim(-5e20, 5e20)
    plt.xlim(-0.05, 0.05)
    plt.savefig(result+'All_Gain{0}.png'.format(suffix))
    #plt.savefig(main_folder+'histograms/{0}.fits'.format(field))


    

   
if combined == True:
    plt.figure(figsize=(10,7))
    plt.subplot2grid((2,2),(0,0))
    x = offset_160
    y = mu_sub_160
    ye = mu_sub_160_sigma

    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit1 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit1)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    
    '''
    std1= np.std(np.delete(y,13))
    # residual and rms only for the blob not including Serpens (index = 13)
    y_new = np.delete(y,13)
    x_new = np.delete(x,13)
    '''
    residual1 = y[small] - (x[small]*fit1[0]+fit1[1])
    rms1 = np.sqrt(np.mean(residual1**2))
    
    #plt.scatter(x[small], y[small], color=c[0], alpha=0.6,ec='w', label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4), round(rms1,4)))
    #plt.scatter(x[big], y[big], marker='*',fmt='o', color=c[0])
    plt.errorbar(x[small], y[small], yerr = ye[small],color=c[0], alpha=0.6,fmt='o', label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4), round(rms1,4)))
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o', color=c[0])

    plt.plot(xp, p(xp), '-', color='k')
    plt.title('160')
    plt.legend(loc=2,ncol=1)
    plt.ylabel(r'Offset $\mu$ from histogram')
    #plt.ylim(-5, 60)

    plt.subplot2grid((2,2),(0,1))
    x = offset_250
    y = mu_sub_250
    ye = mu_sub_250_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit2 = np.polyfit(x[small], y[small], 1)

    p = np.poly1d(fit2)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    '''
    std2= np.std(np.delete(y,13))

    y_new = np.delete(y,13)
    x_new = np.delete(x,13)
    '''
    residual2 = y[small] - (x[small]*fit2[0]+fit2[1])
    rms2 = np.sqrt(np.mean(residual2**2))
    
    #print np.sqrt(np.mean((np.mean(np.delete(y,13))-np.delete(y,13))**2))
    #print np.sqrt(np.mean((np.mean(y)-y)**2))
    plt.errorbar(x[small], y[small], yerr = ye[small], alpha=0.6,fmt='o', color=c[1],label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit2[0],4), round(fit2[1],4), round(rms2,4)))
    plt.plot(xp, p(xp), '-', color='k')
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o', color=c[1])
    #plt.ylim(-20.0, 3.0)
    plt.title('250')
    plt.legend(loc=2,ncol=1)



    plt.subplot2grid((2,2),(1,0))
    x = offset_350
    y = mu_sub_350
    ye = mu_sub_350_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit3 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit3)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    '''
    std3= np.std(np.delete(y,13))

    y_new = np.delete(y,13)
    x_new = np.delete(x,13)
    '''
    residual3 = y[small] - (x[small]*fit3[0]+fit3[1])
    rms3 = np.sqrt(np.mean(residual3**2))
    
    plt.errorbar(x[small], y[small], yerr = ye[small], alpha=0.6,fmt='o', color=c[2],label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit3[0],4), round(fit3[1],4), round(rms3,4)))
    plt.plot(xp, p(xp), '-', color='k')
    plt.errorbar(x[big], y[big], yerr = ye[big],marker='*',fmt='o', color=c[2])
    plt.title('350')
    plt.legend(loc=2,ncol=1)
    plt.ylabel(r'Offset $\mu$ from histogram')
    plt.xlabel(r'Zero-point Correction')


    plt.subplot2grid((2,2),(1,1))
    x = offset_500
    y = mu_sub_500
    ye = mu_sub_500_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]

    fit4 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit4)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    '''
    std4= np.std(np.delete(y,13))
    y_new = np.delete(y,13)
    x_new = np.delete(x,13)
    '''
    residual4 = (y[small] - x[small]*fit4[0]+fit4[1])
    rms4 = np.sqrt(np.mean(residual4**2))
    
    plt.errorbar(x[small], y[small], yerr = ye[small], alpha=0.6,fmt='o', color=c[3],label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit4[0],4), round(fit4[1],4), round(rms4,4)))
    plt.plot(xp, p(xp), '-', color='k')
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o',color=c[3])
    plt.title('500')
    plt.legend(loc=2,ncol=1)
    plt.xlabel(r'Zero-point Correction')

    plt.savefig(result+'Offset_vs_zeropoint{0}.png'.format(suffix))


    mean_line_1 = np.mean((x[small]*fit1[0]+fit1[1]))
    mean_line_2 = np.mean((x[small]*fit2[0]+fit2[1]))
    mean_line_3 = np.mean((x[small]*fit3[0]+fit3[1]))
    mean_line_4 = np.mean((x[small]*fit4[0]+fit4[1]))
#fits = [fit1,fit2,fit3,fit4]
#std = [std1,std2,std3,std4]
#rms = [rms1,rms2,rms3,rms4]
if individual == True:

    plt.figure(figsize=(10,7))
    ax = plt.subplot2grid((1,1),(0,0))
    x = offset_160
    y = mu_sub_160
    ye = mu_sub_160_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit1 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit1)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    #std1= np.std(np.delete(y,13))
    # residual and rms only for the blob not including Serpens (index = 13)
    #y_new = np.delete(y,13)
    #x_new = np.delete(x,13)
    residual1 = y - (x*fit1[0]+fit1[1])
    rms1 = np.sqrt(np.mean(residual1**2))
    ax.errorbar(x[small], y[small], yerr = ye[small], color=c[0],fmt='o',label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4),round(rms1,5)))
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o',color=c[0])
    for i, txt in enumerate(name):
        ax.annotate(txt, (x[i], y[i]))
    plt.plot(xp, p(xp), '-', color='k')
    plt.title('160')
    plt.legend(loc=2,ncol=1)

    plt.savefig(result+'Offset_vs_zeropoint_160{0}.png'.format(suffix))

    plt.figure(figsize=(10,7))
    ax = plt.subplot2grid((1,1),(0,0))
    x = offset_250
    y = mu_sub_250
    ye = mu_sub_250_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit1 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit1)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    #std1= np.std(np.delete(y,13))
    # residual and rms only for the blob not including Serpens (index = 13)
    #y_new = np.delete(y,13)
    #x_new = np.delete(x,13)
    residual1 = y - (x*fit1[0]+fit1[1])
    rms1 = np.sqrt(np.mean(residual1**2))
    ax.errorbar(x[small], y[small], yerr = ye[small], color=c[1],fmt='o',label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4),round(rms1,5)))
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o',color=c[1])
    for i, txt in enumerate(name):
        ax.annotate(txt, (x[i], y[i]))
    plt.plot(xp, p(xp), '-', color='k')
    plt.title('250')
    plt.legend(loc=2,ncol=1)

    plt.savefig(result+'Offset_vs_zeropoint_250{0}.png'.format(suffix))

    plt.figure(figsize=(10,7))
    ax = plt.subplot2grid((1,1),(0,0))
    x = offset_350
    y = mu_sub_350
    ye = mu_sub_350_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit1 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit1)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    #std1= np.std(np.delete(y,13))
    # residual and rms only for the blob not including Serpens (index = 13)
    #y_new = np.delete(y,13)
    #x_new = np.delete(x,13)
    residual1 = y - (x*fit1[0]+fit1[1])
    rms1 = np.sqrt(np.mean(residual1**2))
    ax.errorbar(x[small], y[small], yerr = ye[small], color=c[2],fmt='o',label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4),round(rms1,5)))
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o',color=c[2])
    for i, txt in enumerate(name):
        ax.annotate(txt, (x[i], y[i]))
    plt.plot(xp, p(xp), '-', color='k')
    plt.title('350')
    plt.legend(loc=2,ncol=1)

    plt.savefig(result+'Offset_vs_zeropoint_350{0}.png'.format(suffix))

    plt.figure(figsize=(10,7))
    ax = plt.subplot2grid((1,1),(0,0))
    x = offset_500
    y = mu_sub_500
    ye = mu_sub_500_sigma
    mask = (np.isnan(x) == False) & (np.isnan(y) == False)
    x = x[mask]
    y = y[mask]
    fit1 = np.polyfit(x[small], y[small], 1)
    p = np.poly1d(fit1)
    xp = np.linspace(np.min(x)-10, np.max(x)+10, 100)
    #std1= np.std(np.delete(y,13))
    # residual and rms only for the blob not including Serpens (index = 13)
    #y_new = np.delete(y,13)
    #x_new = np.delete(x,13)
    residual1 = y - (x*fit1[0]+fit1[1])
    rms1 = np.sqrt(np.mean(residual1**2))
    ax.errorbar(x[small], y[small], yerr = ye[small], color=c[3],fmt='o',label = 'y = {0}x + {1}, scatter = {2}'.format(round(fit1[0],4), round(fit1[1],4),round(rms1,5)))
    plt.errorbar(x[big], y[big], yerr = ye[big], marker='*',fmt='o',color=c[3])
    for i, txt in enumerate(name):
        ax.annotate(txt, (x[i], y[i]))
    plt.plot(xp, p(xp), '-', color='k')
    plt.title('500')
    plt.legend(loc=2,ncol=1)

    plt.savefig(result+'Offset_vs_zeropoint_500{0}.png'.format(suffix))



if full == True:

    plt.figure(figsize = (10,7))
    plt.subplot(4,3,1)
    plt.plot(mu_160[small], mu_250[small], '.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('160')
    plt.subplot(4,3,2)
    plt.plot(mu_160[small], mu_350[small], '.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('160')
    plt.subplot(4,3,3)
    plt.plot(mu_160[small], mu_500[small], '.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('160')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,4)
    plt.plot(mu_250[small], mu_160[small], '.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('250')
    plt.subplot(4,3,5)
    plt.plot(mu_250[small], mu_350[small], '.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('250')
    plt.subplot(4,3,6)
    plt.plot(mu_250[small], mu_500[small], '.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('250')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,7)
    plt.plot(mu_350[small], mu_160[small], '.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('350')
    plt.subplot(4,3,8)
    plt.plot(mu_350[small], mu_250[small], '.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('350')
    plt.subplot(4,3,9)
    plt.plot(mu_350[small], mu_500[small], '.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('350')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,10)
    plt.plot(mu_500[small], mu_160[small], '.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('500')
    plt.subplot(4,3,11)
    plt.plot(mu_500[small], mu_250[small], '.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('500')
    plt.subplot(4,3,12)
    plt.plot(mu_500[small], mu_350[small], '.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('500')


    plt.savefig(result+'Correlation_combined{0}.png'.format(suffix))


    plt.figure(figsize = (10,7))
    plt.subplot(4,3,1)
    plt.errorbar(mu_sub_160[small], mu_sub_250[small], xerr = mu_sub_160_sigma[small], yerr=mu_sub_250_sigma[small],fmt='.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('160')
    plt.subplot(4,3,2)
    plt.errorbar(mu_sub_160[small], mu_sub_350[small], xerr = mu_sub_160_sigma[small], yerr=mu_sub_350_sigma[small],fmt='.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('160')
    plt.subplot(4,3,3)
    plt.errorbar(mu_sub_160[small], mu_sub_500[small], xerr = mu_sub_160_sigma[small], yerr=mu_sub_500_sigma[small],fmt='.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('160')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,4)
    plt.errorbar(mu_sub_250[small], mu_sub_160[small], xerr = mu_sub_250_sigma[small], yerr=mu_sub_160_sigma[small],fmt='.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('250')
    plt.subplot(4,3,5)
    plt.errorbar(mu_sub_250[small], mu_sub_350[small], xerr = mu_sub_250_sigma[small], yerr=mu_sub_350_sigma[small],fmt='.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('250')
    plt.subplot(4,3,6)
    plt.errorbar(mu_sub_250[small], mu_sub_500[small], xerr = mu_sub_250_sigma[small], yerr=mu_sub_500_sigma[small],fmt='.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('250')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,7)
    plt.errorbar(mu_sub_350[small], mu_sub_160[small], xerr = mu_sub_350_sigma[small], yerr=mu_sub_160_sigma[small],fmt='.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('350')
    plt.subplot(4,3,8)
    plt.errorbar(mu_sub_350[small], mu_sub_250[small], xerr = mu_sub_350_sigma[small], yerr=mu_sub_250_sigma[small],fmt='.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('350')
    plt.subplot(4,3,9)
    plt.errorbar(mu_sub_350[small], mu_sub_500[small], xerr = mu_sub_350_sigma[small], yerr=mu_sub_500_sigma[small],fmt='.', color=c[2])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('500')
    plt.xlabel('350')

    #plt.figure(figsize = (10,7))
    plt.subplot(4,3,10)
    plt.errorbar(mu_sub_500[small], mu_sub_160[small], xerr = mu_sub_500_sigma[small], yerr=mu_sub_160_sigma[small],fmt='.', color=c[3])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('160')
    plt.xlabel('500')
    plt.subplot(4,3,11)
    plt.errorbar(mu_sub_500[small], mu_sub_250[small], xerr = mu_sub_500_sigma[small], yerr=mu_sub_250_sigma[small],fmt='.', color=c[0])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('250')
    plt.xlabel('500')
    plt.subplot(4,3,12)
    plt.errorbar(mu_sub_500[small], mu_sub_350[small], xerr = mu_sub_500_sigma[small], yerr=mu_sub_350_sigma[small],fmt='.', color=c[1])
    #plt.ylim(-0.05, 0.05)
    #plt.xlim(-0.05, 0.05)
    plt.ylabel('350')
    plt.xlabel('500')


    plt.savefig(result+'Correlation_subtracted_combined{0}.png'.format(suffix))
    