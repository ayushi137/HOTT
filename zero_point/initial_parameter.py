from __future__ import print_function
####################################################################################################
'''
This script has all the Parameters that can be changed. Read the description below
for more information.

By: Ayushi Singh Last Modified: 16 March 2022
'''
###################################################################################################

# path to the folder where everything is. 

path_to_folder = '/mnt/raid-project/hp/asingh/colden_herschel/colden_project'
#path_to_folder = '/Volumes/Storm/Research_has_all_data/colden_project/'


# Prefix of Planck dust model maps. 
# Example TauLL_HFI_CompMap_ThermalDustModel_0_TAN.fits where TauLL is the Planck map name
tau353GHz = '_HFI_CompMap_ThermalDustModel_0_TAN.fits'
temperature = '_HFI_CompMap_ThermalDustModel_4_TAN.fits'
beta = '_HFI_CompMap_ThermalDustModel_6_TAN.fits'

# Dictionary of beam size of each image (in arcsec)
beamsize = {'160':12.5,
			'250':18.4, 
           	'350':25.2,
            '500':36.7, 
            'planck':300.0 
           }

# Frequency of required pixel size ma
regrid_size = 500


########### The following variable are only need for running offset.py ############################

# if need to apply colour correction. This is higly recommended.
apply_cc = True

# if you want to save the colour correction file. It's not required. But true then it
# will be saved in /cc_map/ folder
savecc = True

# allow to fit and add a plane to Herschel Data
add_plane = False

# cropping the high point based on percentage value defined under "percentage". Recommened to be True.
high_crop = True

# high value cropping percentage value if high_crop is True
percentage = 10.0

# cropping the high gradient point based on percentage value defined under "percentage". Recommened to be True.
gradient_crop = True

# high gradient value cropping percentage value if high_crop is True
grad_percentage = 10.0

# cropping the low point based on percentage value defined under "percentage". Recommened to be False.
low_crop = False

# low value cropping percentage value if low_crop is True
low_percentage = 5.0

# cropping based on some residual cutoff. After the first round of fitting the subrated 
# resudual is calculated. Then based on some cut off value the maps is cleaned up again 
# before the next round of fitting. 
residual_crop = True

# sigma clipping needed if the residual_crop is True
sigma_clipping = 3

# Cleaning up the edges so that the clouds are visible and there is not too much white space. 
# The values for edge cutoff are stored in boundary.py under the purpose = edge. 
edge_crop = True

# this value just make the image look "pretty" (for science purpose they do not matter)
# only work with when generating plot for presentaion purpose 
res_lim =0.1

# This make sure that all images are generated with the correlation plot.
images = True

# This allose us to save images and need to be true if you want to save the final result
# can be turned of when just testing the cloud. 
# saved in /herschel_data/<Field>
savefits = True

# saved in /offset_images/
saveimage = True

# allow to write information about offset in a .csv file with today's date in the file name. 
# saved in /offset_values/
write_offset = True 


