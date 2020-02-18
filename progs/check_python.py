

not_found = ''
try:
	import scipy
except:
	not_found=not_found+' scipy'

try:
	import numpy
except:
	not_found=not_found+' numpy'

try:
	import matplotlib
except:
	not_found=not_found+' matplotlib'

if not_found != '':
	print(not_found)

	
