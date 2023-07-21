####################################################################################################
'''
This script make all the needed folders

By: Ayushi Singh
Last Modified: 30 December 2019
'''
###################################################################################################

# required module
import os
from initial_parameter import path_to_folder

### for the offset pipeline

os.system('mkdir '+path_to_folder+'/herschel_archive')

os.system('mkdir '+path_to_folder+'/herschel_data')

os.system('mkdir '+path_to_folder+'/Planck')

os.system('mkdir '+path_to_folder+'/Models')

os.system('mkdir '+path_to_folder+'/grid')

os.system('mkdir '+path_to_folder+'/coverage')

os.system('mkdir '+path_to_folder+'/offset_maps')

os.system('mkdir '+path_to_folder+'/offset_values')

os.system('mkdir '+path_to_folder+'/offset_images')

os.system('mkdir '+path_to_folder+'/offset_mask')

os.system('mkdir '+path_to_folder+'/cc_maps')

### for the SED pipeline

os.system('mkdir '+path_to_folder+'/tau_temp')

os.system('mkdir '+path_to_folder+'/histograms')

os.system('mkdir '+path_to_folder+'/mosaic')
