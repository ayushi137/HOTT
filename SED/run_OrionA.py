from __future__ import print_function
import os, time, datetime

names = [
'OrionA_N',
'OrionA_C',
'OrionA_S']

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




