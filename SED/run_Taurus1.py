from __future__ import print_function
import os, time, datetime

names = [
'TauFill',
'TauL1489',
'TauL1517',
'TauL1521',
'TauL1539',
'TauL1544',
'TauL1551']

date = str((datetime.date.today()))

print ('Date suffix:', date)
print ()
print (names)

for i in range(len(names)):
	start_timer = time.time()
	print (names[i])
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	os.system('python creating_tau_and_temp.py {0}'.format( names[i]))
	
	print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

	end_timer = time.time()
	print (names[i],' took ', str(datetime.timedelta(seconds=end_timer-start_timer)))
	print ()




