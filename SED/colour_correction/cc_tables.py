# color correction
'''
temp_table = np.array([10,15,20,30,40,50,75,100])
	beta_table = np.array([-2,-1,1,2])
correction = np.array([[1.495,1.126,1.035,0.993,0.987,0.987,0.991,0.994],
					   [1.319,1.045,0.988,0.974,0.981,0.988,1.001,1.008],
					   [1.083,0.956,0.957,0.997,1.029,1.051,1.082,1.098],
					   [1.009,1.943,0.971,1.037,1.083,1.113,1.155,1.175]])
'''

from astropy.io import fits
import numpy as np 

a= fits.open('SCalPhotColorCorrK_extended_v5_1449750571140.fits')
table1 = a[1].data
table2 = a[2].data # beta = 0.00
table3 = a[3].data # 0.5
table4 = a[4].data # 1
table5 = a[5].data # 1.25
table6 = a[6].data # 1.5 
table7 = a[7].data # 1.75
table8 = a[8].data # 2.00
table9 = a[9].data # 2.5
table10 = a[10].data # 3.0

Temperature = np.arange(3,300,1)
Beta = np.array([0.,0.5,1.0,1.25,1.5,1.75,2.0,2.5,3.0])

PSW = np.zeros((len(Beta),len(Temperature)))
PMW = np.zeros((len(Beta),len(Temperature)))
PLW = np.zeros((len(Beta),len(Temperature)))

for i in range(2,11,1):
	table = a[i].data
	k = i-2
	for j in range(len(Temperature)):
		PSW[k][j] = table[j][1]
		PMW[k][j] = table[j][2]
		PLW[k][j] = table[j][3] 

temp_table = np.array([10,15,20,30,40,50,75,100])
beta_table = np.array([-2,-1,1,2])
correction_160 = np.array([[1.495,1.126,1.035,0.993,0.987,0.987,0.991,0.994],
					   [1.319,1.045,0.988,0.974,0.981,0.988,1.001,1.008],
					   [1.083,0.956,0.957,0.997,1.029,1.051,1.082,1.098],
					   [1.009,0.943,0.971,1.037,1.083,1.113,1.155,1.175]])


np.savetxt('temperature.txt', Temperature)	
np.savetxt('beta.txt', Beta)	
np.savetxt('250_cc.txt', PSW)		
np.savetxt('350_cc.txt', PMW)		
np.savetxt('500_cc.txt', PLW)
np.savetxt('temperature_160.txt', temp_table)	
np.savetxt('beta_160.txt', beta_table)	
np.savetxt('160_cc.txt', correction_160)			



