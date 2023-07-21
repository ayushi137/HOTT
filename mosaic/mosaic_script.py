from __future__ import print_function
import os, sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rcParams
from astropy import constants as const
import glob
from image_plotting import image_ra_dec
from image_plotting import fit_image_plotting
import required_functions as rf
from turbo_colormap import turbo_cmap
import time, datetime
import colour_correction as cc
import pandas as pd

# alternative to tight layout
rcParams.update({'figure.autolayout': True})
plt.ion()

muh2=2.8                                        # [amu]                 
kappa0 = 10.
h,k,c,mh = const.h.cgs , const.k_B.cgs , const.c.cgs, const.m_p.cgs*1.00794
h,k,c,mh = h.value, k.value, c.value, mh.value  # [erg.s], [erg/K], [cm/s], [g]


start = time.time()

#main_folder = '/Volumes/Storm/Research_has_all_data/colden_project/'
main_folder = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project/'



suffix = '_1_v5' # suffix for the SED folder
suffix_folder = '_v5' # suffix for the mosaic folder
date = '2022-08-21'
suffix2 = '_0_v5' # suffix for the offset values

offset_file = main_folder+'/offset_values/Offset_and_mu_{1}{0}.csv'.format(suffix2,date)

initial = False
frequency = False
SED = True
apply_cc = True

def getoffset(fieldname):
    df = pd.read_csv(offset_file,delimiter=',')
    y = df.index[df['Field']==fieldname]
    print ('The index for ', fieldname, ' is ', y[0])
    offset_values = np.array([df.iloc[y[0]]['Offset 160'],df.iloc[y[0]]['Offset 250'],df.iloc[y[0]]['Offset 350'],df.iloc[y[0]]['Offset 500']])
    return offset_values 

######################### applying colour correction to planck model #########################
def applying_cc(planck, wavelength, temp_map, beta_map):
    print ('Start the colour correction script.')
    # get the temperature and beta from planck maps
    temperature = temp_map
    beta = beta_map

    # get colour correction value from the colour_correction.py
    colour_corrections = np.zeros([np.shape(temperature)[0],np.shape(temperature)[1]])
    for i in range(np.shape(beta)[0]):
        for j in range(np.shape(beta)[1]):
            colour_corrections[i][j] = cc.colourcorrection(temperature[i][j], beta[i][j], wavelength = int(wavelength))
            
    # correct planck 
    planck = planck/colour_corrections
    print ('Finish colour correcting Planck data.\n')
    return planck


field = sys.argv[1]#'Coalsack'

try:
    if field == 'Aquila_Serpens':
        planck_name = 'SerAquL'
        names = ['Aquila','Ser_Main','Serpens','Aquila_W']

    if field == 'Cepheus':
        planck_name = 'CEPHEUSt'
        names = ['CepL1157','CepL1172','CepL1228','CepL1241','CepL1251']

    if field == 'Chameleon':
        planck_name = 'Cha'
        names = ['Cha_I','Cha_II','Cha_III']

    if field == 'Coalsack':
        planck_name = 'Coal'
        names = ['Csack_Glob1','Csack_Glob2','Coalsack']

    if field == 'Corona_Australis':
        planck_name = 'CorAusL'
        names = ['CrA_N','CrA_S']

    if field == 'IC5146':
        planck_name = 'IC5146'
        names = ['IC5146']

    if field == 'Lupus':
        planck_name = 'Lup'
        names = ['Lupus_I','Lupus_III','Lupus_IV-SP2','Lupus_IV-SP1']

    if field == 'Musca':
        planck_name = 'Mus'
        names = ['Musca']

    if field == 'Ophiuchus':
        planck_name = 'OphLL'
        names = ['OphL1688','OphL1712','North_Streamer']

    if field == 'Orion':
        planck_name = 'OriABLL'
        names =['OrionA_N','OrionA_C','OrionA_S','OrionB_N','OrionB_NN','OrionB_S']

    if field == 'Perseus':
        planck_name = 'PerLL'
        names = ['Perseus_E', 'Perseus_W']

    if field == 'Pipe':
        planck_name = 'Pipe'
        names = ['B59','B68','Pipe_C','Pipe_E','Pipe_fill_1','Pipe_fill_2','Pipe_fill_3']

    if field == 'Taurus':
        planck_name = 'TauLL'
        names = ['TauFill','TauL1489','TauL1517','TauL1521','TauL1539','TauL1544','TauL1551','TauS1','TauS2','TauS3','TauT3','TauTMC','TauTMC_E']

    if field == 'W3_HiGal':
        planck_name = 'W3'
        names = ['W3_2']

    if field == 'W3':
        planck_name = 'W3'
        names = ['W3']

    if field == 'Spider':  
        planck_name = 'SPIDERt'
        names = ['Spider']

    if field == 'California':
       planck_name = 'California'
       names = 'California_E'

    if field == 'Draco':  
        planck_name = 'DRACO6t'
        names = ['Draco']    

    if field == 'CygnusX':  
        planck_name = 'CygX'
        names = ['CygX-N','CygX-S']

    if field == 'M16M17':  
        planck_name = 'M16M17'
        names = ['M16','M17']

    if field == 'M16M17_HiGal':  
        planck_name = 'M16M17'
        names = ['M16_2','M17_2']             

    if field == 'MonOB1':  
        planck_name = 'MonOB1'
        names = ['MonOB1']
    
    if field == 'MonR2':  
        planck_name = 'MonR2'
        names = ['MonR2']

    if field == 'NGC2264':  
        planck_name = 'N2264'
        names = ['NGC2264']
        
    if field == 'NGC7538':  
        planck_name = 'N7538'
        names = ['NGC7538']

    if field == 'Rosette':  
        planck_name = 'Rosette'
        names = ['Rosette']

    if field == 'W48':  
        planck_name = 'W48'
        names = ['W48']
        
    if field == 'Vela':  
        planck_name = 'VelaC'
        names = ['Vela']        

    print ('Field:', field, '\nPlanck:', planck_name, '\nNames:',names)

except:
    sys.exit('This field name does not exist. Please check the spelling or modify the script.')





os.system('mkdir {0}mosaic/{1}'.format(main_folder, field+suffix_folder))

#print (planck_name, names)

# Get all the file names
print (main_folder+'mosaic/{1}/{0}_HFI_Model_160.resamp.fits'.format(planck_name, field+suffix_folder))
print (offset_file)


def mygrad(x): # return the absolute value of the gradient 
    diffx,diffy = np.gradient(x)
    return np.sqrt(diffx**2+diffy**2)

planck_orig_tau = main_folder+'Planck_nomaskedyesthreshold/{0}_HFI_CompMap_ThermalDustModel_0_TAN.fits'.format(planck_name)
planck_orig_temp = main_folder+'Planck_nomaskedyesthreshold/{0}_HFI_CompMap_ThermalDustModel_4_TAN.fits'.format(planck_name)
planck_orig_beta = main_folder+'Planck_nomaskedyesthreshold/{0}_HFI_CompMap_ThermalDustModel_6_TAN.fits'.format(planck_name)

if frequency:
    planck_160 = main_folder+'Models_nomaskedyesthreshold/{0}_HFI_model_160.fits'.format(planck_name)
    planck_160_resamp = main_folder+'mosaic/{1}/{0}_HFI_model_160.resamp.fits'.format(planck_name, field+suffix_folder)
    planck_250 = main_folder+'Models_nomaskedyesthreshold/{0}_HFI_model_250.fits'.format(planck_name)
    planck_250_resamp = main_folder+'mosaic/{1}/{0}_HFI_model_250.resamp.fits'.format(planck_name, field+suffix_folder)
    planck_350 = main_folder+'Models_nomaskedyesthreshold/{0}_HFI_model_350.fits'.format(planck_name)
    planck_350_resamp = main_folder+'mosaic/{1}/{0}_HFI_model_350.resamp.fits'.format(planck_name, field+suffix_folder)
    planck_500 = main_folder+'Models_nomaskedyesthreshold/{0}_HFI_model_500.fits'.format(planck_name)
    planck_500_resamp = main_folder+'mosaic/{1}/{0}_HFI_model_500.resamp.fits'.format(planck_name, field+suffix_folder)


planck_tau_n = main_folder+'mosaic/{1}/{0}_tau.fits'.format(planck_name, field+suffix_folder)
planck_temp_n = main_folder+'mosaic/{1}/{0}_temp.fits'.format(planck_name, field+suffix_folder)
planck_colden_n = main_folder+'mosaic/{1}/{0}_colden.fits'.format(planck_name, field+suffix_folder)
planck_beta_n = main_folder+'mosaic/{1}/{0}_beta.fits'.format(planck_name, field+suffix_folder)

planck_tau_resamp = main_folder+'mosaic/{1}/{0}_tau.resamp.fits'.format(planck_name, field+suffix_folder)
planck_temp_resamp = main_folder+'mosaic/{1}/{0}_temp.resamp.fits'.format(planck_name, field+suffix_folder)
planck_colden_resamp = main_folder+'mosaic/{1}/{0}_colden.resamp.fits'.format(planck_name, field+suffix_folder)
planck_beta_resamp = main_folder+'mosaic/{1}/{0}_beta.resamp.fits'.format(planck_name, field+suffix_folder)


template_name = main_folder+'mosaic/{0}/{0}_template.fits'.format(field+suffix_folder)

if frequency:
    data_160_mosaic = main_folder+'mosaic/{0}/{0}_160_mosaic.fits'.format(field+suffix_folder)
    data_160_combined = main_folder+'mosaic/{0}/{0}_160_combined.fits'.format(field+suffix_folder)
    data_250_mosaic = main_folder+'mosaic/{0}/{0}_250_mosaic.fits'.format(field+suffix_folder)
    data_250_combined = main_folder+'mosaic/{0}/{0}_250_combined.fits'.format(field+suffix_folder)
    data_350_mosaic = main_folder+'mosaic/{0}/{0}_350_mosaic.fits'.format(field+suffix_folder)
    data_350_combined = main_folder+'mosaic/{0}/{0}_350_combined.fits'.format(field+suffix_folder)
    data_500_mosaic = main_folder+'mosaic/{0}/{0}_500_mosaic.fits'.format(field+suffix_folder)
    data_500_combined = main_folder+'mosaic/{0}/{0}_500_combined.fits'.format(field+suffix_folder)

grid_mosaic = main_folder+'mosaic/{0}/{0}_grid.fits'.format(field+suffix_folder)
grid_outline = main_folder+'mosaic/{0}/{0}_outline.fits'.format(field+suffix_folder)

if SED:
    tau_mosaic = main_folder+'mosaic/{0}/{0}_tau_mosaic.fits'.format(field+suffix_folder)
    tau_combined = main_folder+'mosaic/{0}/{0}_tau_combined.fits'.format(field+suffix_folder)
    temp_mosaic = main_folder+'mosaic/{0}/{0}_temp_mosaic.fits'.format(field+suffix_folder)
    temp_combined = main_folder+'mosaic/{0}/{0}_temp_combined.fits'.format(field+suffix_folder)
    colden_mosaic = main_folder+'mosaic/{0}/{0}_colden_mosaic.fits'.format(field+suffix_folder)
    colden_combined = main_folder+'mosaic/{0}/{0}_colden_combined.fits'.format(field+suffix_folder)


grid_names = []
data_160_names = [] 
data_250_names = []
data_350_names = []
data_500_names = []
tau_names = []
temp_names = []
colden_names = []
offsets = []

for i in names:
    print (i)
    if frequency:
        a = glob.glob(main_folder+'herschel_data/{0}/{0}*160.offset.conv.resamp.fits'.format(i))
        b = glob.glob(main_folder+'herschel_data/{0}/{0}*250.offset.conv.resamp.fits'.format(i))
        c = glob.glob(main_folder+'herschel_data/{0}/{0}*350.offset.conv.resamp.fits'.format(i))
        d = glob.glob(main_folder+'herschel_data/{0}/{0}*500.offset.fits'.format(i))

        if initial == False:
            offset_values = getoffset(i)
            print ('offset_values:', offset_values)
            offsets.append(offset_values)
        else:
            print ('This is the initial run.') 
            offsets.append(np.array([0.0,0.0,0.0,0.0]))


    g = glob.glob(main_folder+'grid/{0}_500.fits'.format(i))
    if SED:
        x = sorted(glob.glob(main_folder+'tau_temp/{0}{1}/{0}_hott_tau1thz_orig.fits'.format(i, suffix)))
        y = sorted(glob.glob(main_folder+'tau_temp/{0}{1}/{0}_hott_temperature_orig.fits'.format(i, suffix)))
        z = sorted(glob.glob(main_folder+'tau_temp/{0}{1}/{0}_hott_columndensity_orig.fits'.format(i, suffix)))
    print (x)    
    if frequency:
        data_160_names.append(a[0])
        data_250_names.append(b[0])
        data_350_names.append(c[0])
        data_500_names.append(d[0])
    if SED:    
        tau_names.append(x[0])
        temp_names.append(y[0])
        colden_names.append(z[0])
        
    grid_names.append(g[0])
    

temp , head= fits.getdata(planck_orig_temp,  header = True)
beta = fits.getdata(planck_orig_beta)

if SED:
    planck_tau= fits.getdata(planck_orig_tau)
    Tau = planck_tau/((353.e9/1000.e9)**beta)
    CD = planck_tau/(muh2 * mh * (kappa0 / 100.) * ((353.e9/1000.e9)**beta))

    os.system("rm {0}".format(planck_tau_n))
    fits.writeto(planck_tau_n, Tau , head)

    os.system("rm {0}".format(planck_colden_n))
    fits.writeto(planck_colden_n, CD , head)

os.system("rm {0}".format(planck_temp_n))
fits.writeto(planck_temp_n, temp , head)

os.system("rm {0}".format(planck_beta_n))
fits.writeto(planck_beta_n, beta , head)

# get the template map
planck_temp, head = fits.getdata(planck_orig_temp, header = True)

factor = abs(head['CDELT1']*3600/14)
factor_int = int(factor)

head['NAXIS1'] =(int(np.shape(planck_temp)[1]*factor)+1)
head['NAXIS2'] =(int(np.shape(planck_temp)[0]*factor)+1)

template = np.zeros([head['NAXIS2'],head['NAXIS1']])

head['CRPIX1'] =(int(np.shape(planck_temp)[1]*factor/2))
head['CRPIX2'] =(int(np.shape(planck_temp)[0]*factor/2))
head['FIELD'] = (field, 'Name of the field')
head['CDELT1'] = (head['CDELT1']/factor)
head['CDELT2'] = (head['CDELT2']/factor)

try: 
    head.remove(keyword='BUNIT')
    head.remove(keyword='BTYPE')
    head.remove(keyword='TELESCOP')
except:
    print

os.system("rm {0}".format(template_name))
fits.writeto(template_name, template , head)


# regrid planck
if frequency:
    os.system("rm {0}".format(planck_160_resamp))
    rf.regrid(planck_160,template_name, resultimage = planck_160_resamp , header = None)

    os.system("rm {0}".format(planck_250_resamp))
    rf.regrid(planck_250,template_name, resultimage = planck_250_resamp , header = None)

    os.system("rm {0}".format(planck_350_resamp))
    rf.regrid(planck_350,template_name, resultimage = planck_350_resamp , header = None)

    os.system("rm {0}".format(planck_500_resamp))
    rf.regrid(planck_500,template_name, resultimage = planck_500_resamp , header = None)

if SED:
    os.system("rm {0}".format(planck_tau_resamp))
    rf.regrid(planck_tau_n,template_name, resultimage = planck_tau_resamp , header = None)

    os.system("rm {0}".format(planck_colden_resamp))
    rf.regrid(planck_colden_n,template_name, resultimage = planck_colden_resamp , header = None)

os.system("rm {0}".format(planck_temp_resamp))
rf.regrid(planck_temp_n,template_name, resultimage = planck_temp_resamp , header = None)

os.system("rm {0}".format(planck_beta_resamp))
rf.regrid(planck_beta_n,template_name, resultimage = planck_beta_resamp , header = None)

# load template 

template, head = fits.getdata(template_name, header = True)

x = np.shape(template)[0]
y = np.shape(template)[1]

file_num = len(names)


data_160 = np.zeros([file_num, x, y])
data_250 = np.zeros([file_num, x, y])
data_350 = np.zeros([file_num, x, y])
data_500 = np.zeros([file_num, x, y])
result_tau = np.zeros([file_num, x, y])
result_temp = np.zeros([file_num, x, y])
result_colden = np.zeros([file_num, x, y])
grid = np.zeros([file_num, x, y])

# load data files and the grid
print ('loading the data files plus grid')
for i in range(file_num):
    print (i+1, 'of', file_num)

    if frequency:

        offset_val = offsets[i]
        print ('offset values:', offset_val)

        print (data_160_names[i])
        os.system("rm {0}".format(data_160_names[i][:-5]+'.resamp.fits'))
        rf.regrid(data_160_names[i],template_name, resultimage = 'same' , header = None)
        data_160[i]= fits.getdata(data_160_names[i][:-5]+'.resamp.fits') + offset_val[0]

        print (data_250_names[i])
        os.system("rm {0}".format(data_250_names[i][:-5]+'.resamp.fits'))
        rf.regrid(data_250_names[i],template_name, resultimage = 'same' , header = None)
        data_250[i]= fits.getdata(data_250_names[i][:-5]+'.resamp.fits') + offset_val[1]

        print (data_350_names[i])
        os.system("rm {0}".format(data_350_names[i][:-5]+'.resamp.fits'))
        rf.regrid(data_350_names[i],template_name, resultimage = 'same' , header = None)
        data_350[i]= fits.getdata(data_350_names[i][:-5]+'.resamp.fits') + offset_val[2]

        print (data_500_names[i])
        os.system("rm {0}".format(data_500_names[i][:-5]+'.resamp.fits'))
        rf.regrid(data_500_names[i],template_name, resultimage = 'same' , header = None)
        data_500[i]= fits.getdata(data_500_names[i][:-5]+'.resamp.fits') + offset_val[3]
        
    if SED:    
        print (tau_names[i])
        os.system("rm {0}".format(tau_names[i][:-5]+'.resamp.fits'))
        rf.regrid(tau_names[i],template_name, resultimage = 'same' , header = None)
        result_tau[i]= fits.getdata(tau_names[i][:-5]+'.resamp.fits')

        print (temp_names[i])
        os.system("rm {0}".format(temp_names[i][:-5]+'.resamp.fits'))
        rf.regrid(temp_names[i],template_name, resultimage = 'same' , header = None)
        result_temp[i]= fits.getdata(temp_names[i][:-5]+'.resamp.fits')
        
        print (colden_names[i])
        os.system("rm {0}".format(colden_names[i][:-5]+'.resamp.fits'))
        rf.regrid(colden_names[i],template_name, resultimage = 'same' , header = None)
        result_colden[i]= fits.getdata(colden_names[i][:-5]+'.resamp.fits')

    print (grid_names[i])
    os.system("rm {0}".format(grid_names[i][:-5]+'.resamp.fits'))
    rf.regrid(grid_names[i],template_name, resultimage = 'same' , header = None)
    grid[i] = fits.getdata(grid_names[i][:-5]+'.resamp.fits')


grid[grid <= 0.9] = float('nan')
grid[grid >= 1.1] = float('nan')
grid[grid <1.1] = 1.
grid[grid >0.9] = 1.

if frequency:
    data_160_masked = data_160/grid
    data_250_masked = data_250/grid
    data_350_masked = data_350/grid
    data_500_masked = data_500/grid
if SED:
    tau_masked = result_tau/grid
    temp_masked = result_temp/grid
    colden_masked = result_colden/grid
print ('done masking')


if frequency:
# create the herschel mosaic
    data_160_m = np.nanmean(data_160_masked, axis=0)
    os.system("rm {0}".format(data_160_mosaic))
    fits.writeto(data_160_mosaic, data_160_m , head)

    data_250_m = np.nanmean(data_250_masked, axis=0)
    os.system("rm {0}".format(data_250_mosaic))
    fits.writeto(data_250_mosaic, data_250_m , head)

    data_350_m = np.nanmean(data_350_masked, axis=0)
    os.system("rm {0}".format(data_350_mosaic))
    fits.writeto(data_350_mosaic, data_350_m , head)

    data_500_m = np.nanmean(data_500_masked, axis=0)
    os.system("rm {0}".format(data_500_mosaic))
    fits.writeto(data_500_mosaic, data_500_m , head)

grid_m = np.nanmean(grid, axis=0)
os.system("rm {0}".format(grid_mosaic))
fits.writeto(grid_mosaic, grid_m , head)

if SED:
    tau_m = np.nanmean(tau_masked, axis=0)
    os.system("rm {0}".format(tau_mosaic))
    fits.writeto(tau_mosaic, tau_m , head)

    temp_m = np.nanmean(temp_masked, axis=0)
    os.system("rm {0}".format(temp_mosaic))
    fits.writeto(temp_mosaic, temp_m , head)

    colden_m = np.nanmean(colden_masked, axis=0)
    os.system("rm {0}".format(colden_mosaic))
    fits.writeto(colden_mosaic, colden_m , head)


cloudmask_region = (grid_m>0)
gradient_abs = mygrad(cloudmask_region.astype('float64'))
edge_region= (gradient_abs>0)
edge = (edge_region.astype('float64'))
#edge[edge <=0] = float('nan')

os.system('rm {0}'.format(grid_outline))
fits.writeto(grid_outline, edge, head)

print ('mosaic maps created ')

x = np.shape(template)[0]
y = np.shape(template)[1]

# combine herschel mosaic with planck



planck_temp = fits.getdata(planck_temp_resamp)
planck_beta = fits.getdata(planck_beta_resamp)

if SED:
    planck_colden = fits.getdata(planck_colden_resamp)
    planck_tau = fits.getdata(planck_tau_resamp)

if frequency:
    planck_160 = fits.getdata(planck_160_resamp)
    if apply_cc == True:
        print ('Start the colour correction script.')
        planck_160 = applying_cc(planck_160,160, planck_temp, planck_beta)    
        print ('Finish colour correcting Planck data.\n')

    planck_250 = fits.getdata(planck_250_resamp)
    if apply_cc == True:
        print ('Start the colour correction script.')
        planck_250 = applying_cc(planck_250,250, planck_temp, planck_beta)    
        print ('Finish colour correcting Planck data.\n')

    planck_350 = fits.getdata(planck_350_resamp)
    if apply_cc == True:
        print ('Start the colour correction script.')
        planck_350 = applying_cc(planck_350,350, planck_temp, planck_beta)    
        print ('Finish colour correcting Planck data.\n')

    planck_500 = fits.getdata(planck_500_resamp)
    if apply_cc == True:
        print ('Start the colour correction script.')
        planck_500 = applying_cc(planck_500,500, planck_temp, planck_beta)    
        print ('Finish colour correcting Planck data.\n')


print ('combining data')
if frequency:
    combined_160 = np.zeros([x, y])
    combined_250 = np.zeros([x, y])
    combined_350 = np.zeros([x, y])
    combined_500 = np.zeros([x, y])

if SED:
    combined_tau = np.zeros([x, y])
    combined_temp = np.zeros([x, y])
    combined_colden = np.zeros([x, y])

print (np.shape(planck_temp))
for i in range(x):
    print (i+1, 'of total', x)
    for j in range(y):
        if np.isnan(grid_m[i][j]) == True:
            if frequency:
                combined_160[i][j] = planck_160[i][j]
                combined_250[i][j] = planck_250[i][j]
                combined_350[i][j] = planck_350[i][j]
                combined_500[i][j] = planck_500[i][j]

            if SED:
                combined_tau[i][j] = planck_tau[i][j]
                combined_temp[i][j] = planck_temp[i][j]
                combined_colden[i][j] = planck_colden[i][j]
                
        else:
            if frequency:
                combined_160[i][j] = data_160_m[i][j]
                combined_250[i][j] = data_250_m[i][j]
                combined_350[i][j] = data_350_m[i][j]
                combined_500[i][j] = data_500_m[i][j]

            if SED:
                combined_tau[i][j] = tau_m[i][j]
                combined_temp[i][j] = temp_m[i][j]
                combined_colden[i][j] = colden_m[i][j]

if frequency:
    os.system("rm {0}".format(data_160_combined))
    fits.writeto(data_160_combined, combined_160 , head)  

    os.system("rm {0}".format(data_250_combined))
    fits.writeto(data_250_combined, combined_250 , head)  

    os.system("rm {0}".format(data_350_combined))
    fits.writeto(data_350_combined, combined_350 , head)  

    os.system("rm {0}".format(data_500_combined))
    fits.writeto(data_500_combined, combined_500 , head) 

if SED:
    os.system("rm {0}".format(tau_combined))
    fits.writeto(tau_combined, combined_tau , head) 

    os.system("rm {0}".format(temp_combined))
    fits.writeto(temp_combined, combined_temp , head) 

    os.system("rm {0}".format(colden_combined))
    fits.writeto(colden_combined, combined_colden , head) 

print ("done")

print ("\n************ This field name was:", field, "*********************")
totalruntime = str(datetime.timedelta(seconds=time.time()-start))
print ('\nRun Time: ', str(datetime.timedelta(seconds=time.time()-start)))


try:
    os.system('echo Mosaic finished for the field {0} in {1}. | mailx -r asingh@cita.utoronto.ca -s "Code Finished" ayushi.singh@mail.utoronto.ca'.format(field, totalruntime))
except:
    pass
