##############################################################################
'''
This script has a function that store all the offset values generated from the
histogram analysis of the Initial run of the data.  

This script requires following functions and files:
initial_parameter.py

The programs consists of following functions:
- getoffset 			: Store all the offset values 

By: Ayushi Singh
Last Modified: 18 March 2022 
'''
##############################################################################
import numpy as np
from initial_parameter import offset_file
import pandas as pd 

def getOffset(fieldname):
	df = pd.read_csv(offset_file,delimiter=',')
	y = df.index[df['Field']==fieldname]
	print ('The index for ', fieldname, ' is ', y[0])
	offset_values = np.array([df.iloc[y[0]]['Offset 160'],df.iloc[y[0]]['Offset 250'],df.iloc[y[0]]['Offset 350'],df.iloc[y[0]]['Offset 500']])
	return offset_values
