from __future__ import print_function

# This script is there to run all the fields at once. 
# Only use this if you know what each script does and the data has been set up 
# using rename.py and make_folders.py

import os, time, datetime
from initial_parameter import path_to_folder


#### Gould Belt regions (Polaris is missing)

names = [
'Aquila',
'Ser_Main',
'Serpens',
'Aquila_W',
'CepL1157',
'CepL1172',
'CepL1228',
'CepL1241',
'CepL1251',
'Cha_I',
'Cha_II',
'Cha_III',
'Csack_Glob1',
'Csack_Glob2',
'Coalsack',
'CrA_N',
'CrA_S',
'IC5146',
'Lupus_I',
'Lupus_III',
'Lupus_IV-SP2',
'Lupus_IV-SP1',
'Musca',
'OphL1688',
'OphL1712',
'North_Streamer',
'OrionA_N',
'OrionA_C',
'OrionA_S',
'OrionB_N',
'OrionB_NN',
'OrionB_S',
'Perseus_E',
'Perseus_W',
'Pipe_C',
'Pipe_E',
'B59',
'B68',
'Pipe_fill_1',
'Pipe_fill_2',
'Pipe_fill_3',
'TauFill',
'TauL1489',
'TauL1517',
'TauL1521',
'TauL1539',
'TauL1544',
'TauL1551',
'TauS1',
'TauS2',
'TauS3',
'TauT3',
'TauTMC',
'TauTMC_E',
'Spider',
'Draco',
'California_E',
'CygX-N',
'CygX-S',
'M16',
'M17',
'MonOB1',
'MonR2',
'NGC2264',
'NGC7538',
'Rosette',
'W48',
'W3',
'Vela',
'W3_2',
'M16_2',
'M17_2']

planck = [ 
'SerAquL',
'SerAquL',
'SerAquL',
'SerAquL',
'CEPHEUSt',
'CEPHEUSt',
'CEPHEUSt',
'CEPHEUSt',
'CEPHEUSt',
'Cha',
'Cha',
'Cha',
'Coal',
'Coal',
'Coal',
'CorAusL',
'CorAusL',
'IC5146',
'Lup',
'Lup',
'Lup',
'Lup',
'Mus',
'OphLL',
'OphLL',
'OphLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'PerLL',
'PerLL',
'Pipe',
'Pipe',
'Pipe',
'Pipe',
'Pipe',
'Pipe',
'Pipe',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'SPIDERt',
'DDRACO6t',
'California',
'CygX',
'CygX',
'M16M17',
'M16M17',
'MonOB1',
'MonR2',
'N2264',
'N7538',
'Rosette',
'W48',
'W3',
'VelaC',
'W3',
'M16M17',
'M16M17' ] 


'''
# HOBYS regions that are slow maps, Fast maps for M16, M17 and W3, Vela and Spider 
names = [ 'CygX-N',
'CygX-S',
'M16',
'M17',
'MonOB1',
'MonR2',
'NGC2264',
'NGC7538',
'Rosette',
'W48',
'W3',
'Vela',
'W3_2',
'M16_2',
'M17_2',
'Spider']

planck = [ 'CygX',
'CygX',
'M16M17',
'M16M17',
'MonOB1',
'MonR2',
'N2264',
'N7538',
'Rosette',
'W48',
'W3',
'VelaC',
'W3',
'M16M17',
'M16M17',
'Spider'] 


names = ['Aquila',
'Ser_Main',
'Serpens',
'Aquila_W',
'CrA_N',
'CrA_S',
'IC5146',
'OphL1688',
'OphL1712',
'North_Streamer',
'OrionA_N',
'OrionA_C',
'OrionA_S',
'OrionB_N',
'OrionB_NN',
'OrionB_S',
'Perseus_E',
'Perseus_W',
'Spider',
'W3_2',
'TauFill',
'TauL1489',
'TauL1517',
'TauL1521',
'TauL1539',
'TauL1544',
'TauL1551',
'TauS1',
'TauS2',
'TauS3',
'TauT3',
'TauTMC',
'TauTMC_E']

planck = ['SerAquL',
'SerAquL',
'SerAquL',
'SerAquL',
'CorAusL',
'CorAusL',
'IC5146',
'OphLL',
'OphLL',
'OphLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'PerLL',
'PerLL',
'SPIDERt',
'W3',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL']

names = ['California_E',
'California_W']
planck =  ['California',
'California']



names = [
'M16',
'M17',
'MonOB1',
'MonR2',
'NGC2264',
'NGC7538',
'Rosette',
'W48',
'W3',
'Vela',
'W3_2',
'M16_2',
'M17_2']

planck = [ 
'M16M17',
'M16M17',
'MonOB1',
'MonR2',
'N2264',
'N7538',
'Rosette',
'W48',
'W3',
'VelaC',
'W3',
'M16M17',
'M16M17' ] 
'''

names = [
'OphL1688',
'OphL1712',
'North_Streamer',
'OrionA_N',
'OrionA_C',
'OrionA_S',
'OrionB_N',
'OrionB_NN',
'OrionB_S',
'Perseus_E',
'Perseus_W',
'TauFill',
'TauL1489',
'TauL1517',
'TauL1521',
'TauL1539',
'TauL1544',
'TauL1551',
'TauS1',
'TauS2',
'TauS3',
'TauT3',
'TauTMC',
'TauTMC_E',]

planck = [ 
'OphLL',
'OphLL',
'OphLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'OriABLL',
'PerLL', 
'PerLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL',
'TauLL'
] 

names  = ['W3', 'W3_2']
planck = ['W3', 'W3']
print (names, planck)
for i in range(len(names)):
	start_timer = time.time()
	print (names[i], planck[i])
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	
	print ('\nRunning Backward SED Modeling')
	#os.system('python backward_sed_modelling.py {0}'.format(planck[i]))
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	
	print ('\nRunning Unit Change')
	#os.system('python unit_change.py {0}'.format(names[i]))
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	
	print ('\nRunning Creating CSR')
	#os.system('python creating_csr.py {0} {1}'.format(planck[i], names[i]))
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	print ('\nRunning Regrid Planck')
	#os.system('python regrid_planck.py {0} {1}'.format(planck[i], names[i]))
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	
	print ('\nRunning Zero Point Correction')
	os.system('python zero_point_correction.py {0} {1} 160'.format(planck[i], names[i]))
	os.system('python zero_point_correction.py {0} {1} 250'.format(planck[i], names[i]))
	os.system('python zero_point_correction.py {0} {1} 350'.format(planck[i], names[i]))
	os.system('python zero_point_correction.py {0} {1} 500'.format(planck[i], names[i]))
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	
	end_timer = time.time()
	totalruntime = str(datetime.timedelta(seconds=end_timer-start_timer))
	print (names[i],' took ', totalruntime)
	print ()
	
	try:
		os.system('echo zero-point correction for the field {0} in {1}. | mailx -r asingh@cita.utoronto.ca -s "Code Finished" ayushi.singh@mail.utoronto.ca'.format(names[i], totalruntime))
	except:
		pass
	

dateToday = (datetime.date.today()) 
os.system('mkdir {0}/offset_images/offsets-{1}/before'.format(path_to_folder,str(dateToday)))
os.system('mkdir {0}/offset_images/offsets-{1}/after'.format(path_to_folder,str(dateToday)))

os.system('mv {0}/offset_images/offsets-{1}/*before.png {0}/offset_images/offsets-{1}/before/'.format(path_to_folder,str(dateToday)))
os.system('mv {0}/offset_images/offsets-{1}/*after.png {0}/offset_images/offsets-{1}/after/'.format(path_to_folder,str(dateToday)))

