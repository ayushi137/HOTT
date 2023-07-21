from __future__ import print_function
import os, time, datetime

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
'W3',
'W3_2',
'Spider',
'California_E',
'California_W']

'''
names = [
'OrionA_N',
'OrionA_C',
'OrionA_S',
'OrionB_N',
'OrionB_NN',
'OrionB_S',
'Perseus_E',
'Perseus_W',
'TauL1521']
'''


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

names=[
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
'TauTMC_E']

names = ['W3', 'W3_2']
date = str((datetime.date.today()))

print ('Date suffix:', date)
print ()
print ('Total number of files:', len(names))
print (names)
for i in range(len(names)):
	start_timer = time.time()
	print (names[i])
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	os.system('python post_SED_hist_good.py {0} {1}'.format(names[i], date))
	print ()

	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	end_timer = time.time()
	print (names[i],' took ', str(datetime.timedelta(seconds=end_timer-start_timer)))
	print ()

