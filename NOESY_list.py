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


# cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
cwd = '/Users/mclay1/FGFR3_structure/cyana_23/'
calc = cwd + 'CALC.cya'
cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
noa = cwd + 'cycle7.noa'

print(cya_plists)

noalines = open(noa).readlines()
pl = 0
unusedupl = open('final_unused_3.upl','w')
usedupl = open('final_used.upl','w')
Used_upls = []
Used_aupls = []
for noelist in cya_plists:
	pl+=1
	nota,single,amb,notused,peak,incr,dia,ndia  = 0, 0, 0, 0, 0, 0, 0,0
	for x in range(len(noalines)):
		line = noalines[x]
		if noelist in line and 'out of' in noalines[x+1]:
			peak+=1
			peakn = line.strip().split()[1]
			dist = line.strip().split()[-2]
			if '0 out of 0' in noalines[x+1]:
				nota+=1
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				notused+=1
			if '1 out of' in noalines[x+1] and 'diagonal' in line:
				single+= 1
			if '0 out of' not in noalines[x+1] and 'diagonal' not in line:
				single+= 1
				# print(noalines[x+ int(noalines[x+1].split()[3])+2])
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*'] and cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if noalines[x+y].strip()[0] not in ['!','*'] and cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					upl = '%4d %-4s %-4s %4d %-4s %-4s' %(resi1, resn1, atom1, resi2, resn2, atom2)
					upl2 = '%4d %-4s %-4s %4d %-4s %-4s' %(resi2, resn2, atom2, resi1, resn1, atom1)
					if upl not in Used_upls and pshift > 0.9 and maxd < 8.0:
						# print(upl)
						Used_aupls.extend([upl,upl2])
						usedupl.write('%4d %-4s %-4s %4d %-4s %-4s    %s  #peak %s #plist %d #Pshift %0.2f\n' %(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift ))
			if noalines[x+1].strip().split()[0] > '1' and 'diagonal' not in line:
				amb+=1
			if 'increased' in line:
				incr+=1 
			if 'diagonal' in line and '0 out of' not in noalines[x+1]:
				dia+=1
			if 'diagonal' in line and '0 out of' in noalines[x+1]:
				ndia+=1
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if noalines[x+y].strip()[0] not in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					upl = '%4d %-4s %-4s %4d %-4s %-4s' %(resi1, resn1, atom1, resi2, resn2, atom2)
					upl2 = '%4d %-4s %-4s %4d %-4s %-4s' %(resi2, resn2, atom2, resi1, resn1, atom1)
					if upl not in Used_upls and pshift > 0.9 and maxd < 8.0:
						# print(upl)
						Used_upls.extend([upl,upl2])
						unusedupl.write('%4d %-4s %-4s %4d %-4s %-4s    %s  #peak %s #plist %d #Pshift %0.2f\n' %(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift ))
	print(noelist)
	print("Total number of peaks %d" %peak)
	print('Number with single assignment %d'%single)
	print('Number ambigious assignment %d'%amb)
	print('Number assigned but unused %d'%notused)
	print('Number never assinged %d'%nota)
	print('Number of diagonal assignments %d' %dia)
	print('Number of diagonal not assignmened %d' %ndia)
	print('Number of peaks with increase distacne cut off %d' %incr)
	print('total %d' %(nota+notused+single+amb+dia+ndia))
	print()
unusedupl.close()

