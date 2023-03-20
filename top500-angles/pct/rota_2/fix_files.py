import os
import numpy as np
files = ['rota500-arg.data',
'rota500-asn.data',
'rota500-asp.data',
'rota500-gln.data',
'rota500-glu.data',
'rota500-his.data',
'rota500-ile.data',
'rota500-leu.data',
'rota500-lys.data',
'rota500-met.data',
'rota500-phetyr.data',
'rota500-trp.data']
# for i in range(0,360,10):
# 	print(i)
# for i in np.arange(5):
# 	print(i)

# intervlas = [np.arrange(0,15,5),np.arrange(0,10,1)]
# for ran in intervlas:


# for file in files:
# 	print(file)
# 	ndim = open(file).readlines()[1][-2]
# 	print(ndim)
# 	increments = []
# 	for x in range(3,int(ndim)+3,1):
# 		increments.append([0.0, float(open(file).readlines()[x].split()[-3]),float(open(file).readlines()[x].split()[-3])/int(open(file).readlines()[x].split()[-2])])
# 	print(increments)
# count = 5 
# trans_dict = {}
# for i in np.arange(0,360,5):
# 	for x in np.arange(0,360,5):
# 		count+=1
# 		for ni in np.arange(i,i+5,1):
# 			for nx in np.arange(x,x+5,1):
# 				print(ni,nx, count)
				
# print(count)
MQ = {}
count = 6 
trans_dict = {}
for i in np.arange(0,360,8):
	for x in np.arange(0,360,8):
		for y in np.arange(0,360,8):
			count+=1
			for ni in np.arange(i,i+8,1):
				for nx in np.arange(x,x+8,1):
					# for ny in np.arange(x,x+8,1):
					MQ[count] = [ni,nx]
					print(ni,nx, count)
# print(count)