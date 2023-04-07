'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
# 11/10/2022
#
# Updated January 30, 2023 to include peak intensities 
# Updated March 09, 2023 to account for old cycle7.peaks format missing information for multiple assignment crosspeaks
# Updated March 14, 2023 to add header to output list, and short constraints tag
# ------------------------------------------------------------------------------

'''
import os
import sys
import numpy as np
import re
import glob

AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "HIST": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S'}


def analize_noa(cwd, outdir, calc, noa7, Seqdict, violdict, qupldict,upldict,pad,upldict2):
	cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
	prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
	log = glob.glob(os.path.join(cwd + 'log*'))[0]

	plist_dict = {}
	for x in range(len(cya_plists)):
		plist = cya_plists[x].replace('.peaks','')
		plist_dict[plist] = str(x)
		exec("unused{:} = ['### Peaks which CYANA did not use assignment possiblity #unused\\n']".format(str(x)))
		exec("no_assign{:} = ['### Peaks which CYANA found no assignment possibility #no assignmnet\\n']".format(str(x)))
		exec("good{:} = open('{:}{:}.list','w')".format(str(x), outdir, plist))
		exec("good{:}.write('### Peaks which CYANA assigned\\n#{:}  {:^26}  {:^24}  {:^5}  {:^10}  {:^6}  {:^24}\\n')".format(str(x),'Peak #','Frequencies','Connection','UPL', 'Range' ,'Pshift','Comment'))
		exec("questlist{:} = []".format(str(x)))
		exec("usedquestlist{:} = []".format(str(x)))
		exec("peaks{:} = {{}}".format(str(x)))
		exec("intensity{:} = {{}}".format(str(x)))
	for x in range(len(cya_plists)):
		intdict = eval('intensity' + str(x))
		pdict = eval('peaks' + str(x))
		peaklines = open(cwd + cya_plists[x].replace('.peaks','-cycle7.peaks')).readlines()
		for i in range(len(peaklines)):
			line = peaklines[i]
			if line.strip():
				if line.strip()[0] != '#':
					if line[0:7] == '       ':  # account for old format of cycle7.peaks
						line = peaklines[i-1][0:peaklines[i-1].find(' U ')+35] + ' ' + peaklines[i].strip()
						peaklines[i] = line
					pdict[int(line.split()[0])] = [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
					if line.split()[5] == 'U':  # get intensity for 3D NOESY
						intensity = "{:12.6E}".format(float(line.split()[6]))
						if len(line.split()) > 13 and line.split()[13] == '#VC':
							intensity = "{:12.6E}".format(float(line.split()[6]) * float(line.split()[14]))
					elif line.split()[6] == 'U':  # get intensity of psudo 4D NOESY
						intensity = "{:12.6E}".format(float(line.split()[7]))
						if len(line.split()) > 15 and line.split()[15] == '#VC':
							intensity = "{:12.6E}".format(float(line.split()[7]) * float(line.split()[16]))
					if int(line.split()[0]) in intdict.keys():
						intdict[int(line.split()[0])].append(intensity)
					else:
						intdict[int(line.split()[0])] = [intensity]

	Calibration_cns = []
	Swapped = {}
	for line in open(log).readlines():
		if "Calibration constant for peak list" in line:
			newline = line.replace(line.split()[5],cya_plists[int(line.split()[5][0:-1]) -1].replace('.peaks',':'))
			if newline not in Calibration_cns:
				Calibration_cns.append(newline)
		if line.strip() and line.strip().split()[-1] == 'swapped':
			res = line.strip().split()
			group1 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[2])
			group2 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[3])
			Swapped[group1] = group2
			Swapped[group2] = group1
	Assignments= []
	for prot in prots:
		for line in open(cwd + prot.replace('.prot','-final.prot')).readlines():
			if line.strip():
				if '#' != line.strip()[0]:
					resi = line.strip().split()[4]
					atom = line.strip().split()[3]
					assign = Seqdict[resi] + '-' + atom
					if assign not in Assignments and atom[0] in ['H','Q']:
						Assignments.append(assign)
	Assignments = sorted(Assignments, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))

	assigndict,assigndict2 = {}, {}
	from itertools import combinations
	ADpairs = [ '{:}-{:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
	for comb in ADpairs:
		assigndict[comb] = []
		assigndict2[comb] = []
	noalines = open(noa7).readlines()
	for x in range(len(noalines)):
		line = noalines[x]
		if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
			if '0 out of' not in noalines[x+1] and 'diagonal' not in noalines[x]:
				peak = int(noalines[x].strip().split()[1])
				plist = noalines[x].strip().split()[3].replace('.peaks','')
				pdict = eval('peaks' + plist_dict[plist])
				intdict = eval('intensity' + plist_dict[plist])
				questionable = eval('questlist' + plist_dict[plist])
				used =eval('usedquestlist' + plist_dict[plist])
				SUP = noalines[x+1].split()[-1].replace(':','')
				linepad = pad[len(plist):]
				Calconst = float(Calibration_cns[int(plist_dict[plist])].split()[-1])
				# print(plist)
				# print('{:1.3E}'.format(Calconst))
				for y in range(2,int(noalines[x+1].split()[0])+2,1):
					cns = noalines[x+y].strip().split()
					if cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5], cns[6], int(cns[7]), float(cns[10])/100, cns[13]
					if cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[0],cns[1],int(cns[2]), cns[4], cns[5], int(cns[6]), float(cns[9])/100, cns[12]
					group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
					group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
					note = ''
					if 'increased' in line: note = note + ' increased '
					if group1 in Swapped.keys(): group1 = Swapped[group1]; note = note + ' swapped '
					if group2 in Swapped.keys(): group2 = Swapped[group2]; note = note + ' swapped '
					conect = '{:}-{:}'.format(group1,group2)
					if '{:} peak {:} from {:}'.format(conect,peak,plist) in upldict2: 
						note = ' UPL ' + note
					if float(intdict[peak][y-2]) != 0.0: 
						dist = (Calconst/float(intdict[peak][y-2]))**(1/6)
						drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
					if float(intdict[peak][y-2]) == 0.0: 
						print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
					outline = '{:^28} {:^14} {:>9}A  Peak {:4} from {:<}{:}  pshift {:3.2f} {:}\n'.format(conect,intdict[peak][y-2],drange,peak,linepad,plist,pshift,note)
					udist = ' '
					if conect in upldict.keys(): udist = upldict[conect] + 'A'
					if '{:}-{:}'.format(group1,group2)in ADpairs:
						assigndict2['{:}-{:}'.format(group1,group2)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
						assigndict['{:}-{:}'.format(group1,group2)].append(outline)
					if '{:}-{:}'.format(group2,group1)in ADpairs:
						assigndict['{:}-{:}'.format(group2,group1)].append(outline)
						assigndict2['{:}-{:}'.format(group2,group1)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
					if float(SUP) <= 0.6:
						questionable.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   low support\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',float(SUP)))
						used.append(peak)
					if conect in violdict.keys():
						outline = "{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:6.2f}   {:}".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,violdict[conect].replace(' #',''))
						questionable.append(outline)
						used.append(peak)
					if conect in qupldict.keys():
						outline = "{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   {:}".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,qupldict[conect])
						questionable.append(outline)
						used.append(peak)
					if pshift <= 0.60 and peak not in used:
						questionable.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   poor chem shift\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist,drange +'A',pshift))
						used.append(peak)
	print('finished assigned')

	ADpairs2 = []
	for con in ADpairs:
		if len(assigndict[con]) >= 1:
			ADpairs2.append(con)
	for x in range(len(noalines)):
		line = noalines[x]
		if 'Peak' in noalines[x] and '0 out of' in noalines[x+1]:
			peak = int(noalines[x].strip().split()[1])
			plist = noalines[x].strip().split()[3].replace('.peaks','')
			pdict = eval('peaks' + plist_dict[plist])
			intdict = eval('intensity' + plist_dict[plist])
			noassingmnet = eval('no_assign' + plist_dict[plist])
			unused = eval('unused' + plist_dict[plist])
			linepad = pad[len(plist):]
			if '0 out of 0' not in noalines[x+1]:
				nopt = int(noalines[x+1].split()[3])
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if len(cns) > 8:
						if noalines[x+y].strip()[0] in ['!','*']:
							atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
						if noalines[x+y].strip()[0] not in ['!','*']:
							atom1,resn1, resi1, atom2, resn2, resi2,pshift, drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
						group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
						group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
						if group1 in Swapped.keys(): group1 = Swapped[group1]
						if group2 in Swapped.keys(): group2 = Swapped[group2]
						conect = '{:}-{:}'.format(group1,group2)
						if float(intdict[peak][0]) != 0.0: 
							dist = (Calconst/float(intdict[peak][0]))**(1/6)
							drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
						if float(intdict[peak][0]) == 0.0: 
							print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
						outline = '#{:^28} {:^14} {:>9}A  Peak {:4} from {:<}{:}  pshift {:0.2f} unused\n'.format(conect,intdict[peak][0],drange,peak,plist,linepad,pshift)
						if '{:}-{:}'.format(group1,group2)in ADpairs2:
							assigndict['{:}-{:}'.format(group1,group2)].append(outline)
							assigndict2['{:}-{:}'.format(group1,group2)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
						if '{:}-{:}'.format(group2,group1)in ADpairs2:
							assigndict['{:}-{:}'.format(group2,group1)].append(outline)
							assigndict2['{:}-{:}'.format(group2,group1)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
						udsit = ' '
						if conect in upldict.keys(): udist = upldict[conect]+'A'
						if pshift > 0.75 and conect in upldict.keys():
							unused.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   unused {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,noalines[x+nopt+2].strip()[:-1].replace('Violated','Viol').replace('structures ','')))
						if pshift > 0.75 and conect not in upldict.keys():
							unused.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   unused \n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist,drange +'A',pshift))
						if pshift < 0.75 and int(noalines[x+1].split()[3]) == 1:
							noassingmnet.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   unused {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,noalines[x+nopt+2].strip()[:-1].replace('Violated','Viol').replace('structures ','')))
			if '0 out of 0' in noalines[x+1]:
				drange = ''
				if float(intdict[peak][0]) != 0.0: 
					dist = (Calconst/float(intdict[peak][0]))**(1/6)
					drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
				if float(intdict[peak][0]) == 0.0: 
					print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
				pdict = eval('peaks' + plist_dict[plist])
				noassingmnet.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6}   no assignmnet\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2], 'none','', drange+'A','na'))


	assigned = open(outdir + 'Assignment_Summary.txt','w')
	for val in Calibration_cns:
		assigned.write(val.strip()+'\n')
	assigned.write('\n\n')
	for con in ADpairs:
		if len(assigndict[con]) >= 1:
			if con in upldict.keys():
				assigned.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
			if con not in upldict.keys():
				assigned.write('{:} ({:}):\n'.format(con,len(assigndict[con])))
			assigned.writelines(assigndict[con])
			assigned.write('\n')
		if len(assigndict[con]) == 1:
			pline = assigndict[con][0].split()
			peak = int(pline[4])
			plist = pline[6]
			pshift = pline[8]
			pdict = eval('peaks' + plist_dict[plist])
			questout =  eval('questlist' + plist_dict[plist])
			used =eval('usedquestlist' + plist_dict[plist])
			if peak not in used:
				upl = ' '
				note = 'only '
				if con in upldict.keys(): upl = upldict[con] + 'A'
				if con in qupldict.keys(): note = note + qupldict[con]
				questout.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6}   {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],pline[0],upl, pline[2],pline[8],note))
				used.append(peak)

	for x in range(len(noalines)):
		line = noalines[x]
		if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
			if '0 out of' not in noalines[x+1] and 'diagonal' not in noalines[x]:
				peak = int(noalines[x].strip().split()[1])
				plist = noalines[x].strip().split()[3].replace('.peaks','')
				pdict = eval('peaks' + plist_dict[plist])
				intdict = eval('intensity' + plist_dict[plist])
				used =eval('usedquestlist' + plist_dict[plist])
				SUP = noalines[x+1].split()[-1].replace(':','')
				linepad = pad[len(plist):]
				good = eval('good' + plist_dict[plist])
				Calconst = float(Calibration_cns[int(plist_dict[plist])].split()[-1])
				for y in range(2,int(noalines[x+1].split()[0])+2,1):
					cns = noalines[x+y].strip().split()
					if cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5], cns[6], int(cns[7]), float(cns[10])/100, cns[13]
					if cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[0],cns[1],int(cns[2]), cns[4], cns[5], int(cns[6]), float(cns[9])/100, cns[12]
					group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
					group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
					note = ''
					if group1 in Swapped.keys(): group1 = Swapped[group1]; note = note + 'swapped '
					if group2 in Swapped.keys(): group2 = Swapped[group2]; note = note + 'swapped '
					conect = '{:}-{:}'.format(group1,group2)
					if '{:} peak {:} from {:}'.format(conect,peak,plist) in upldict2: 
						note = note + 'UPL '
					if float(intdict[peak][y-2]) != 0.0: 
						dist = (Calconst/float(intdict[peak][y-2]))**(1/6)
						drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
					if float(intdict[peak][y-2]) == 0.0: 
						print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
					udist = ''
					if conect in upldict.keys(): udist = upldict[conect] + 'A'
					if peak not in used:
						good.write("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,note))

	for x in range(len(cya_plists)):
		eval("good{:}.close()".format(str(x)))
		exec("good{:} = open('{:}{:}.list','a')".format(str(x), outdir, cya_plists[x].replace('.peaks','')))
		unused = eval('unused{:}'.format(str(x)))
		noassingmnet = eval('no_assign{:}'.format(str(x)))
		exec("good{:}.write('### Peaks with SUP < 0.6 #low support\\n### Peaks with poor chemical shift matching (Pshift < 0.6) #poor chem shift\\n### Peak associated with UPL > 6.0 #long distance\\n### Peak associated with UPL < 3.0 #short distance\\n### Peak represening only instance of an assignment #only\\n')".format(str(x)))
		questlist = eval('questlist' + plist_dict[cya_plists[x].replace('.peaks','')])
		qlist = sorted(questlist, key = lambda x: float(x.strip().split()[0]))
		eval("good{:}.writelines(qlist)".format(str(x)))
		eval("good{:}.writelines(unused)".format(str(x)))
		eval("good{:}.writelines(noassingmnet)".format(str(x)))
		eval("good{:}.close()".format(str(x)))
	assigned.close()
	print("Finished cycle7.noa analysis")

	return assigndict2


