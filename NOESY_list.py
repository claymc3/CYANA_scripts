'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
# 06/28/2022
# ------------------------------------------------------------------------------

'''
import os
import sys
import numpy as np


cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
calc = cwd + 'CALC.cya'
cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
noa = cwd + 'cycle7.noa'

print(cya_plists)

noalines = open(noa).readlines()

for noelist in cya_plists:
	nacount,uucount,peak = 0,0,0
	amb, single = 0, 0
	incr = 0 
	dia = 0
	for x in range(len(noalines)):
		line = noalines[x]
		if noelist in line and 'out of' in noalines[x+1]:
			peak+=1
			if '0 out of 0' in noalines[x+1]:
				nacount+=1
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				uucount+=1
			if '1 out of' in noalines[x+1]:
				single+= 1
			if noalines[x+1].strip().split()[0] > '1':
				amb+=1
			if 'increased' in line:
				incr+=1 
			if 'diagonal' in line:
				dia+=1
	print(noelist)
	print("Total number of peaks %d" %peak)
	print('Number never assinged %d'%nacount)
	print('Number assigned but unused %d'%uucount)
	print('Number with single assignment %d'%single)
	print('Number ambigious assignment %d'%amb)
	print('Number of diagonal assignments %d' %dia)
	print('Number of peaks with increase distacne cut off %d' %incr)
	print('total %d' %(nacount+uucount+single+amb))
	print()