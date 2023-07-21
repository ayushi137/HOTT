from __future__ import print_function
import os, time, datetime

names = [
'B59',
'B68',
'Pipe_C',
'Pipe_E',
'Pipe_fill_1',
'Pipe_fill_2',
'Pipe_fill_3'
]

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




