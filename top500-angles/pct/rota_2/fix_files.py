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


files = ['rota500-arg.data',
'rota500-lys.data']
# 4
# [[0.0, 360.0, 10.0], [0.0, 360.0, 10.0], [0.0, 360.0, 10.0], [0.0, 360.0, 10.0]]
for file in files:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 7
	outfile.writelines(open(file).readlines()[0:7])
	# newlines = [open(file).readlines()[0]]
	lines = open(file).readlines()
	for x in np.arange(0,360,10):
		for y in np.arange(0,360,10):
			xylist = 0.0
			for z in np.arange(0,360,10):
				for n in np.arange(0,360,10):
					count+=1
					xylist = xylist + float(lines[count].split()[-1])
			for nx in np.arange(x,x+10,1):
				for ny in np.arange(y,y+10,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,xylist)
					outfile.write(outline)
		# print(count)
	outfile.close()

files = ['rota500-gln.data',
'rota500-met.data']
for file in ['rota500-gln.data']:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 6
	outfile.writelines(open(file).readlines()[0:6])
	# newlines = [open(file).readlines()[0]]
	lines = open(file).readlines()
	for x in np.arange(0,360,8):
		for y in np.arange(0,360,8):
			xylist = 0.0
			for z in np.arange(0,360,8):
				print(z)
				count+=1
				print('count {:}'.format(count))
				print(lines[count])
				xylist = xylist + float(lines[count].split()[-1])
			for nx in np.arange(x,x+8,1):
				for ny in np.arange(y,y+8,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,xylist/26.3477149248961)
					outfile.write(outline)
		# print(count)
	outfile.close()
for file in ['rota500-met.data']:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 6
	outfile.writelines(open(file).readlines()[0:6])
	# newlines = [open(file).readlines()[0]]
	lines = open(file).readlines()
	for x in np.arange(0,360,8):
		for y in np.arange(0,360,8):
			xylist = 0.0
			for z in np.arange(0,360,8):
				print(z)
				count+=1
				print('count {:}'.format(count))
				print(lines[count])
				xylist = xylist + float(lines[count].split()[-1])
			for nx in np.arange(x,x+8,1):
				for ny in np.arange(y,y+8,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,xylist/18.0087045570916)
					outfile.write(outline)
		# print(count)
	outfile.close()


files = ['rota500-phetyr.data',
'rota500-asp.data']
# 2
# [[0.0, 360.0, 5.0], [0.0, 180.0, 5.0]]

for file in files:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 5 
	outfile.writelines(open(file).readlines()[0:6])
	newlines = []
	lines = open(file).readlines()
	for x in np.arange(0,360,5):
		for y in np.arange(0,180,5):
			count+=1
			for nx in np.arange(x,x+5,1):
				for ny in np.arange(y,y+5,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,lines[count].split()[-1],nx,ny+180,lines[count].split()[-1])
					outfile.write(outline)
	outfile.close()

files = ['rota500-asn.data',
'rota500-his.data',
'rota500-trp.data',
'rota500-ile.data',
'rota500-leu.data']

for file in files:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 5 
	outfile.writelines(open(file).readlines()[0:6])
	lines = open(file).readlines()
	for x in np.arange(0,360,5):
		for y in np.arange(0,360,5):
			count+=1
			for nx in np.arange(x,x+5,1):
				for ny in np.arange(y,y+5,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,lines[count].split()[-1])
					outfile.write(outline)
	outfile.close()
# 2
# [[0.0, 360.0, 5.0], [0.0, 360.0, 5.0]]

# rota500-glu.data
# 3
# [[0.0, 360.0, 8.0], [0.0, 360.0, 8.0], [0.0, 180.0, 7.826086956521739]]
for file in ['rota500-glu.data']:
	outfile = open(file.replace('.data','_chi1chi2.data'),'w')
	count = 6
	outfile.writelines(open(file).readlines()[0:6])
	lines = open(file).readlines()
	for x in np.arange(0,360,8):
		for y in np.arange(0,360,8):
			xylist = 0.0
			for z in range(23):
				count+=1
				xylist = xylist + float(lines[count].split()[-1])
			for nx in np.arange(x,x+8,1):
				for ny in np.arange(y,y+8,1):
					outline = '{:<4.1f} {:<4.1f}  {:}\n'.format(nx,ny,xylist)
					outfile.write(outline)
	outfile.close()

