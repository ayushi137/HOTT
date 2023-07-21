import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
plt.ion()


def image_ra_dec(name, save_image = None, 
                 colour_bar_label = r'$N_{H_{2}}$ log$[cm^{-2}]$', 
                 min_limit = 20.45, max_limit = 22.5, style = 'inferno', 
                 length = 10.5, height = 6.5, log_map = True):
	
	'''
    ##### THINGS YOU NEED TO CHANGE 
    name = 'FILE.fits'	# path and name for the input file
    min_limit = 20.45	# min limit for the colour bar # play with these numbers 
    max_limit = 22.5	# max limit for the colour bar # play with these numbers 
    style = 'inferno'	# image colour theme
    colour_bar_label = '$N_{H_{2}}$ log$[cm^{-2}]$' # label for the colour bar
    save_image = 'file.png' # name and path for the resulting file 
    length = 10.5 		# length of the image # play with these number 
    height = 6.5 		# heigth of the image # play with these number 
    log_map = True          # make a log10 image
    ######
    '''

	import numpy as np
	import matplotlib.pyplot as plt
	from astropy.wcs import WCS
	from astropy.io import fits
	from astropy.utils.data import get_pkg_data_filename

	plt.figure(figsize=(length,height))

	filename = get_pkg_data_filename(name)

	hdu = fits.open(filename)[0]
	wcs = WCS(hdu.header)

	if log_map == True:
	    data_image = np.log10(hdu.data)
	else:
	    data_image = hdu.data
	ax = plt.subplot(projection=wcs)
	im = ax.imshow(data_image, vmin=min_limit, vmax=max_limit, origin='lower', cmap=style)

	ax = plt.gca()
	lon = ax.coords[0]
	lat = ax.coords[1]

	lon.set_major_formatter('hh:mm')
	lat.set_major_formatter('dd:mm')
	lon.set_separator(('h', "m", 's'))
	lon.set_axislabel('Right Ascension (J2000)')
	lat.set_axislabel('Declination (J2000)')
	cbar = plt.colorbar(im)
	cbar.set_label(colour_bar_label,fontsize=12,rotation=90,labelpad=15)
	if save_image == None:
		plt.savefig(name[:-5]+'.png', transparent=True)
	else:
		plt.savefig(save_image, transparent=True)

	return ()


def fit_image_plotting(name, save_image = None, 
                 colour_bar_label = r'$N_{H_{2}}$ log$[cm^{-2}]$', 
                 VMIN = 20.45, VMAX = 22.5, CMAP = 'inferno', 
                 length = 9, height = 6, stretch = 'log', transparent= False):
	import aplpy
	fig = plt.figure(figsize=(length, height))
	f1  = aplpy.FITSFigure(name,figure=fig,subplot=(1,1,1))

	f1.show_colorscale(cmap=CMAP,vmin=VMIN,vmax=VMAX, stretch=stretch)

	# tick coordinates
	f1.tick_labels.set_xformat("hh:mm")
	f1.tick_labels.set_yformat("dd")
	# axis labels
	f1.axis_labels.set_font(size=20)
	# tick labels
	f1.tick_labels.show()
	f1.tick_labels.set_font(size=20)
	# colour bar
	f1.add_colorbar()
	f1.colorbar.set_axis_label_text(colour_bar_label)
	f1.colorbar.set_axis_label_font(size=20)

	f1.add_scalebar(1.0,color="white",corner="top left")
	f1.scalebar.set_label("1 degrees")
	f1.scalebar.set_font(size=20)
	plt.subplots_adjust(top = 1,bottom=0,right=1,left=0,hspace=0,wspace=0)
	fig.canvas.draw()
	f1.save(name.split(".fits")[0]+".pdf")

	return()
