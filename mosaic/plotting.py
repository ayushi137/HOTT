import os,sys
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

#%matplotlib inline
#%matplotlib notebook
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=UserWarning) 
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

main_folder = '/Volumes/Storm/Research_has_all_data/colden_project/'
#main_folder = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project/'

suffix_folder = '_01'#'_planck_wo_plane'
field = sys.argv[1]#'Aquila_Serpens'
frequency =True
SED = True
values = {'Aquila_Serpens':{'length' : 9, 'height' : 6, 'vmax_160' : 700, 'vmax_250' : 500, 'vmax_350' : 350, 'vmax_500' : 200, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.003, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Cepheus':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Chameleon':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 21, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 8e19, 'vmax_cd' : 2e22},
          'Corona_Australis':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 8e19, 'vmax_cd' : 2e22},
          'Coalsack':{'length' : 9, 'height' : 6, 'vmax_160' : 700, 'vmax_250' : 500, 'vmax_350' : 350, 'vmax_500' : 200, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.003, 'vmin_cd' : 2e20, 'vmax_cd' : 2e22},
          'IC5146':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 25, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Lupus':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Musca':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 21, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 1e22},
          'Ophiuchus':{'length' : 9, 'height' : 6, 'vmax_160' : 200, 'vmax_250' : 170, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 27, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Orion':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 25, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Perseus':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Pipe':{'length' : 9, 'height' : 6, 'vmax_160' : 200, 'vmax_250' : 170, 'vmax_350' : 100, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'Taurus':{'length' : 9, 'height' : 6, 'vmax_160' : 150, 'vmax_250' : 120, 'vmax_350' : 80, 'vmax_500' : 30, 
                    'vmin_temp' : 10, 'vmax_temp' : 20, 'vmin_tau' : 0, 'vmax_tau' : 0.001, 'vmin_cd' : 1e20, 'vmax_cd' : 2e22},
          'W3':{'length' : 9, 'height' : 6, 'vmax_160' : 700, 'vmax_250' : 500, 'vmax_350' : 350, 'vmax_500' : 200, 
                    'vmin_temp' : 10, 'vmax_temp' : 27, 'vmin_tau' : 0, 'vmax_tau' : 0.002, 'vmin_cd' : 3e20, 'vmax_cd' : 2e22},
          'W3_HOBYS':{'length' : 9, 'height' : 6, 'vmax_160' : 700, 'vmax_250' : 500, 'vmax_350' : 350, 'vmax_500' : 200, 
                    'vmin_temp' : 10, 'vmax_temp' : 27, 'vmin_tau' : 0, 'vmax_tau' : 0.002, 'vmin_cd' : 3e20, 'vmax_cd' : 2e22},
          'Spider':{'length' : 9, 'height' : 6, 'vmax_160' : 30, 'vmax_250' : 20, 'vmax_350' : 10, 'vmax_500' : 5, 
                    'vmin_temp' : 15, 'vmax_temp' : 23, 'vmin_tau' : 0, 'vmax_tau' : 0.0001, 'vmin_cd' : 1e19, 'vmax_cd' : 4e20},
         }

length = values[field]['length']
height = values[field]['height']
vmax_160 = values[field]['vmax_160']
vmax_250 = values[field]['vmax_250']
vmax_350 = values[field]['vmax_350']
vmax_500 = values[field]['vmax_500']
vmin_temp = values[field]['vmin_temp']
vmax_temp = values[field]['vmax_temp']
vmin_tau = values[field]['vmin_tau']
vmax_tau = values[field]['vmax_tau']
vmin_cd = values[field]['vmin_cd']
vmax_cd = values[field]['vmax_cd']



def fit_image_plotting(name, save_image = None, 
                       label = r'$N_{H_{2}}$ log$[cm^{-2}]$', 
                       VMIN = 20.45, VMAX = 22.5, CMAP = 'inferno', 
                       length = 12, height = 9, stretch = 'log', 
                       transparent= False, scalebar = True, scalebar_label = True,
                      contour= None, levels = 1, outline_color = 'white'):
    import aplpy
    fig = plt.figure(figsize=(length, height))
    f1  = aplpy.FITSFigure(name,figure=fig)
    
    f1.show_colorscale(cmap=CMAP,vmin=VMIN,vmax=VMAX, stretch=stretch) 
    
    if contour != None:
        f1.show_contour(contour,colors=outline_color, levels=levels)

    # tick coordinates
    f1.tick_labels.set_xformat("hh:mm")
    f1.tick_labels.set_yformat("dd")
    # axis labels
    
    f1.axis_labels.set_font(size=15)
    # tick labels
    f1.tick_labels.show()
    f1.tick_labels.set_font(size=15)
    # colour bar
    f1.add_colorbar()
    f1.colorbar.set_axis_label_text(label)
    f1.colorbar.set_axis_label_font(size=20)
    
    if scalebar:
        f1.add_scalebar(1.0,color="white",corner="top left")
        if scalebar_label:
            f1.scalebar.set_label("1 degree")
            f1.scalebar.set_font(size=20)
        plt.subplots_adjust(top = 1,bottom=0,right=1,left=0,hspace=0,wspace=0)
        fig.canvas.draw()
    if save_image == None:    
        f1.save(name.split(".fits")[0]+".png", transparent = transparent)
    else:
        f1.save(save_image, transparent = transparent)
    
    return()
    

field = field+suffix_folder
if frequency:    
  data_160_mosaic = main_folder+'mosaic/{0}/{0}_160_mosaic.fits'.format(field)
  data_160_combined = main_folder+'mosaic/{0}/{0}_160_combined.fits'.format(field)
  data_250_mosaic = main_folder+'mosaic/{0}/{0}_250_mosaic.fits'.format(field)
  data_250_combined = main_folder+'mosaic/{0}/{0}_250_combined.fits'.format(field)
  data_350_mosaic = main_folder+'mosaic/{0}/{0}_350_mosaic.fits'.format(field)
  data_350_combined = main_folder+'mosaic/{0}/{0}_350_combined.fits'.format(field)
  data_500_mosaic = main_folder+'mosaic/{0}/{0}_500_mosaic.fits'.format(field)
  data_500_combined = main_folder+'mosaic/{0}/{0}_500_combined.fits'.format(field)

grid_mosaic = main_folder+'mosaic/{0}/{0}_grid.fits'.format(field)
grid_outline = main_folder+'mosaic/{0}/{0}_outline.fits'.format(field)

if SED:
  tau_mosaic = main_folder+'mosaic/{0}/{0}_tau_mosaic.fits'.format(field)
  tau_combined = main_folder+'mosaic/{0}/{0}_tau_combined.fits'.format(field)
  temp_mosaic = main_folder+'mosaic/{0}/{0}_temp_mosaic.fits'.format(field)
  temp_combined = main_folder+'mosaic/{0}/{0}_temp_combined.fits'.format(field)
  colden_mosaic = main_folder+'mosaic/{0}/{0}_colden_mosaic.fits'.format(field)
  colden_combined = main_folder+'mosaic/{0}/{0}_colden_combined.fits'.format(field)

if frequency:  
  name = data_160_combined
  fit_image_plotting(name, label = r"Intensity at 160$\mu m$ [MJy/sr]", 
               VMIN = 0, VMAX = vmax_160,length = length, height = height, stretch = 'linear')
  print ('Finished 160 Frequency map')

  name = data_250_combined
  fit_image_plotting(name, label = r"Intensity at 250$\mu m$ [MJy/sr]", 
               VMIN = 0, VMAX = vmax_250,length = length, height = height, stretch = 'linear')            
  print ('Finished 250 Frequency map')

  name = data_350_combined
  fit_image_plotting(name, label = r"Intensity at 350$\mu m$ [MJy/sr]", 
               VMIN = 0, VMAX = vmax_250,length = length, height = height, stretch = 'linear')
  print ('Finished 350 Frequency map')

  name = data_500_combined
  fit_image_plotting(name, label = r"Intensity at 500$\mu m$ [MJy/sr]", 
               VMIN = 0, VMAX = vmax_500,length = length, height = height, stretch = 'linear')
  print ('Finished 500 Frequency map')

if SED:
  name = tau_combined
  fit_image_plotting(name, label = r"Optical Depth at 1THz", 
               VMIN = vmin_tau, VMAX = vmax_tau,length = length, height = height, stretch = 'linear')
  print ('Finished tau map')

  name = temp_combined
  fit_image_plotting(name, label = r"Dust Temperature [K]", 
                     VMIN = vmin_temp, VMAX = vmax_temp, 
                     length = length, height = height,stretch = 'linear', CMAP = 'gray')
  print ('Finished temp map')

  name = colden_combined
  fit_image_plotting(name, label = r"$N_{\rm{H}_{2}}$ log$[cm^{-2}]$", 
               VMIN = vmin_cd, VMAX =vmax_cd,length = length, height = height, stretch = 'log', CMAP = 'inferno')
  print ('Finished colden map')

  name = colden_combined
  fit_image_plotting(name, label = r"$N_{\rm{H}_{2}}$ log$[cm^{-2}]$", save_image=name.split(".fits")[0]+"_outlined.png",
               VMIN = vmin_cd, VMAX =vmax_cd,length = length, height = height, stretch = 'log', CMAP = 'inferno', contour = grid_outline)    
  print ('Finished colden with outline map')

  name = temp_combined
  fit_image_plotting(name, label = r"Dust Temperature [K]", save_image=name.split(".fits")[0]+"_outlined.png",
                     VMIN = vmin_temp, VMAX = vmax_temp, 
                     length = length, height = height,stretch = 'linear', CMAP = 'gray', contour = grid_outline, outline_color = 'purple')
  print ('Finished temp with ouline map')

  name = tau_combined
  fit_image_plotting(name, label = r"Optical Depth at 1THz", save_image=name.split(".fits")[0]+"_outlined.png", 
               VMIN = vmin_tau, VMAX = vmax_tau,length = length, height = height, stretch = 'linear',contour = grid_outline)
  print ('Finished tau with outline map')

