### Mary Clay
import os
import sys
import glob
import pandas as pd
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB',
'CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB',
'ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB',
'GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG',
'PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2',
'GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA',
'HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1',
'ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEQD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1',
'LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ',
'LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUQQD':'CD2,CD1',
'METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHG2':'CG','METHG3':'CG','METQG':'CG','METQE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE',
'ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD2':'ND2',
'PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD',
'GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNQE2':'NE2',
'ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2',
'SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG',
'THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2',
'VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALQQG':'CG1,CG2',
'TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2',
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH'}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',
if len(sys.argv)==1:
	print('''

Usage: 
	cyanapra [pdb] [upl]

Required Input:

	PDB			PDB to be used typically the final.pdb or pdb after CNS refinment
				If this is not located in current directory provide path
					CNS/refinePDB/r12_cya.pdb

	upl				What upl file would you like to use? final.upl cycle?.upl
					Which ever upl you specify determines the overveiw file used. 

OutPut:

	name_dist.cxc
	name_dist.pml
	peak_list(s).pb
	name_poor.pb
	name_long.pb
	name_short.pb
	name_viol.pb
''')
	exit()

cwd = os.getcwd() + '/'
outdir = cwd + 'post_cyana_ana/'
print(outdir)
in_pdb = sys.argv[1]
fupl = sys.argv[2]
pdbname = in_pdb.split('.')[0]
fovw = fupl.replace('.upl','.ovw')
calc = cwd + 'CALC.cya'
outname = fupl.split('.')[0]

print(os.path.exists(outdir))
if not os.path.exists(outdir):
	os.makedirs(outdir)

checkcons = open(outdir + outname + '_summary.txt','w')
checkcons.write('                         #peaks    upl Violations Assigned Ambiguous Unassigned\n')
summary = pd.DataFrame(columns=['#peaks', 'upl', 'Violations', 'Assigned', 'Ambiguous', 'Unassigned' ])
cya_plists = [line.strip().replace('.peaks','-cycle7.peaks') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and 'constraints' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
print(upls)
lols = [con for con in manualongcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualongcons if 'aco' in con]
plistdict = {}
for x in range(len(cya_plists)):
	exec('upl' + str(x+1) + ' = []')
	exec('pb' + str(x+1) + ' = []')
	plist = cya_plists[x]
	plistdict[str(x+1)] = plist.replace('-cycle7.peaks','')
	upl = [line.strip() for line in open(fupl).readlines() if line.strip() and 'plist '+ str(x+1) in line]
	viol = [line.strip() for line in open(fovw).readlines() if line.strip() and 'list '+ str(x+1) in line]
	na,sa,aa = 0, 0, 0 
	cyplines = [line for line in open(os.getcwd() + '/' + plist).readlines() if line.strip() and line[0] != "#"]
	for line in cyplines:
		if 'e 0     0     0     0' in line: na+=1
		if len(line[0:8].strip()) !=0 and 'VC' in line: aa+=1
		if 'e 0     0     0     0' not in line and 'VC' not in line: sa+=1
	checkcons.write('%-25s%6d %6d %10d %8d %9d %10d\n' %(plist.replace('-cycle7.peaks',''),na+aa+sa,len(upl),len(viol),sa,aa,na))
	summary.loc[plist.replace('-cycle7.peaks',''),'#peaks'] = na+aa+sa
	summary.loc[plist.replace('-cycle7.peaks',''),'upl'] = len(upl)
	summary.loc[plist.replace('-cycle7.peaks',''),'Violations'] = len(viol)
	summary.loc[plist.replace('-cycle7.peaks',''),'Unassigned'] = na
	summary.loc[plist.replace('-cycle7.peaks',''),'Ambiguous'] = aa
	summary.loc[plist.replace('-cycle7.peaks',''),'Assigned'] = sa
print(summary)
checkcons.write('\n\n')
colors = ['white','raspberry','gold','forest','marine','purple','orange','cyan','pink','deepteal','gray']
colors2 = ['white','mediumvioletred','orange','forest','royalblue','purple','chocolate','cyan','pink','deepteal','gray']

outpml = open(outdir + fupl.replace('.upl','_pra.pml'),'w')
outpml.write('load '+ cwd + in_pdb+'\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('color gray60, all\n')
outcmx = open(outdir + fupl.replace('.upl','_pra.cxc'),'w')
outcmx.write('open '+ cwd + in_pdb+'\n')
pdbname = in_pdb.replace('.pdb','')


badpbout = open(outdir + outname + '_poor_cons.pb','w')
badpbout.write("; halfbond = false\n; color = darkred\n; radius = 0.2\n; dashes = 0\n")
longpbout = open(outdir + outname + '_long_cons.pb','w')
longpbout.write("; halfbond = false\n; color = aquamarine\n; radius = 0.2\n; dashes = 0\n")
shortpbout = open(outdir + outname + '_short_cons.pb','w')
shortpbout.write("; halfbond = false\n; color = light coral\n; radius = 0.2\n; dashes = 0\n")
pviolpbout = open(outdir + outname + '_viol_peaks_cons.pb','w')
pviolpbout.write("; halfbond = false\n; color = brown\n; radius = 0.2\n; dashes = 0\n")
uviolpbout = open(outdir + outname + '_viol_upl_cons.pb','w')
uviolpbout.write("; halfbond = false\n; color = brown\n; radius = 0.2\n; dashes = 0\n")
'''cns[0] = resi1 , cns[1] = resn1, cns[2]= atom1, cns[3] = resi2, cns[4]=resn2, cns[5] = atom2, cns[6] = dist'''
i = 1
finalupl = []
# name1 = line[10:15].strip()
# resn1 = line[16:20]
# resi1 = line[20:24].strip()
# name2 = line[27:32].strip()
# resn2 = line[33:37]
# resi2 = line[37:41].strip()
Upperdict, Lowerdict,= {}, {}
pmlpviols,pmluviols = [],[]
Filtered = []
checkcons.write('### Violated Distance Constraints from %s \n' %(str(fovw)))
v = 1
for line in open(fovw).readlines():
	if line[4:9] == 'Upper' or line[4:9] == 'Lower':
		v+=1
		atom1 = line[10:15].strip()
		if line[16:20].strip()+line[10:15].strip() in replacements.keys():
			atom1 = replacements[line[16:20].strip()+line[10:15].strip()]
		atom2 = line[27:32].strip()
		if line[33:37].strip()+line[27:32].strip() in replacements.keys():
			atom2 = replacements[line[33:37].strip()+line[27:32].strip()]
		if 'peak' not in line:
			cons = '%4s %4s %-3s   %4s %4s %-3s   %6.2f\n' % (line[20:24].strip(),line[16:20],line[10:15].strip(),line[37:41].strip(),line[33:37],line[27:32].strip(),float(line[44:48]))
			cons2 = '%4s %4s %-3s   %4s %4s %-3s   %6.2f  # %s %s\n' % (line[20:24].strip(),line[16:20],line[10:15].strip(),line[37:41].strip(),line[33:37],line[27:32].strip(),float(line[44:48]), line[50:52], line[62:66])
			if line[4:9] == 'Upper':
				Upperdict[cons] = cons2
			if line[4:9] == 'lower':
				Lowerdict[cons] = cons2
			if line[50:52].strip() >= '10':
				uviolpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(line[20:24].strip(), atom1, line[37:41].strip(),atom2, 'brown'))
				dist = 'distance uplviol%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(v), pdbname, line[20:24].strip(), atom1, pdbname, line[37:41].strip(), atom2)
				tcolor = 'color brown, upl viol ' + str(v) + '\n'
				pmluviols.append('uplviol' + str(v))
		if 'peak' in line and 'QQ' not in line:
			pviolpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(line[20:24].strip(), atom1, line[37:41].strip(),atom2, 'brown'))
			dist = 'distance peakviol %s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(v), pdbname, line[20:24].strip(), atom1, pdbname, line[37:41].strip(), atom2)
			tcolor = 'color brown, peakviol ' + str(v) + '\n'
			pmlpviols.append('peakviol' + str(v))
			for line2 in open(fupl).readlines():
				cns = line2.split()
				if cns[8] == line[90:].split()[1] and cns[10] == line[90:].split()[3] and cns[2] == line[10:15].strip() and cns[5] == line[27:32].strip():
					checkcons.write(line2.replace('\n',' #Violated ' + line[50:88]+ '\n'))
					Filtered.append(line2)
			outpml.write(dist)
			outpml.write(tcolor)
checkcons.write('\n\n')
finalupl,poorcons,poorcons2, show, shortcons, shortcons2, longcons,longcons2 = [],[],[],[],[],[],[],[]

for line in open(fupl).readlines():
	if 'QQ' in line.split()[2] or 'QQ' in line.split()[5]:
		pass 
	else:
		cns = line.split()
		upllist = eval('upl' + cns[10])
		pblist = eval('pb' + cns[10])
		atom1 = cns[2]
		atom2 = cns[5]
		if cns[1]+cns[2] in replacements.keys():
			atom1 = atom1.replace(cns[2], replacements[cns[1]+cns[2]])
		if cns[4]+cns[5] in replacements.keys():
			atom2 = atom2.replace(cns[5], replacements[cns[4]+cns[5]])
		if cns[1]+cns[5] not in replacements.keys():
			atom1 = atom1
		if cns[4]+cns[5] not in replacements.keys():
			atom2=atom2
		i+=1
		if len(line.split()) < 12:
			upllist.append('UPL' + str(i))
			dist = 'distance UPL%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
			tcolor = 'color ' + colors[int(cns[10])]  + ', UPL' + str(i) + '\n'
			pblist.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, colors2[int(cns[10])]))
		if len(line.split()) >= 12:
			if float(cns[12]) > 0.5:
				upllist.append('UPL' + str(i))
				dist = 'distance UPL%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
				tcolor = 'color ' + colors[int(cns[10])]  + ', UPL' + str(i) + '\n'
				pblist.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, colors2[int(cns[10])]))
				finalupl.append(line)
			if float(cns[12]) < 0.5:
				badpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'darkred'))
				poorcons2.append(line)
				Filtered.append(line)
				dist = 'distance poor%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
				tcolor = 'color red, poor' + str(i) + '\n'
				poorcons.append('poor'+ str(i))
				# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
				# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
			if float(cns[6]) >= 6.0:
				if line not in Filtered:
					Filtered.append(line)
					longpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'aquamarine'))
					longcons2.append(line)
					dist = 'distance long%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
					tcolor = 'color cyan, long' + str(i) + '\n'
					longcons.append('long' + str(i))
					# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
					# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
			if float(cns[6]) <= 3.00:
				if abs(int(cns[0])- int(cns[3])) > 1:
					# if atom1 != 'H' or atom2 != 'H':
						if line not in Filtered:
							Filtered.append(line)
							shortpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'light coral'))
							shortcons2.append(line)
							dist = 'distance short%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
							tcolor = 'color orange, short' + str(i) + '\n'
							shortcons.append('short' + str(i))
							# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
							# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
		outpml.write(dist)
		outpml.write(tcolor)
badpbout.close()
longpbout.close()
shortpbout.close()
pviolpbout.close()
uviolpbout.close()

for uplfile in upls:
	newlines = []
	fin = open(uplfile,'r')
	for line in fin.readlines():
		if line in Upperdict.keys():
			newlines.append(Upperdict[line])
		if line not in Upperdict.keys():
			newlines.append(line)
	fout = open(uplfile,'w')
	fout.writelines(newlines)
fin.close()
fout.close()
u = 1
for uplfile in upls:
	fin = open(uplfile,'r')
	outpb = open(outdir + uplfile.replace('.upl','_cons.pb'),'w')
	pmlgroup = 'group %s, ' %(uplfile.replace('.upl',''))
	outpb.write("; halfbond = false\n; color = pink\n; radius = 0.1\n; dashes = 10\n")
	for line in fin.readlines():
		cns = line.split()
		if "#" not in cns[0]:
			u+=1
			atom1 = cns[2]
			atom2 = cns[5]
			if cns[1]+cns[2] in replacements.keys():
				atom1 = atom1.replace(cns[2], replacements[cns[1]+cns[2]])
			if cns[4]+cns[5] in replacements.keys():
				atom2 = atom2.replace(cns[5], replacements[cns[4]+cns[5]])
			if cns[1]+cns[5] not in replacements.keys():
				atom1 = atom1
			if cns[4]+cns[5] not in replacements.keys():
				atom2=atom2
			if ',' in atom1 and ',' not in atom2:
				outpb.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1.split(',')[0], cns[3],atom2, 'pink'))
				outpb.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1.split(',')[1], cns[3],atom2, 'pink'))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1.split(',')[0], pdbname, cns[3], atom2))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1.split(',')[1], pdbname, cns[3], atom2))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
			if ',' in atom2 and ',' not in atom1:
				outpb.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2.split(',')[0], 'pink'))
				outpb.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2.split(',')[1], 'pink'))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2.split(',')[0]))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2.split(',')[2]))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
			if ',' not in atom1 and ',' not in atom2:
				outpb.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'pink'))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
	outpml.write(pmlgroup + '\n')

for lolfile in lols:
	newlines = []
	fin = open(lolfile,'r')
	for line in fin.readlines():
		if line in Lowerdict.keys():
			newlines.append(Lowerdict[line])
		if line not in Lowerdict.keys():
			newlines.append(line)
	fout = open(lolfile,'w')
	fout.writelines(newlines)
fin.close()
fout.close()


filtered_upl = open(fupl.replace('.upl','4cns.upl'),'w')
for upl in finalupl:
	filtered_upl.write(upl)
filtered_upl.close()
# print(poorcons2)
for x in range(len(cya_plists)):
	plist = cya_plists[x].replace('-cycle7.peaks','')
	outname = fupl.split('.')[0]
	outcmx.write('open ' + outdir + outname + '_'+ plist + '.pb\n')
	outcmx.write('color #%s %s\n' %(str(x+2),colors2[x+1]))
	pbout = open(outdir + outname + '_'+ plist + '.pb','w')
	pbout.write("; halfbond = false\n")
	pbout.write("; color = " + colors2[x+1] + '\n')
	pbout.write("; radius = 0.1\n")
	pbout.write("; dashes = 0\n")
	upllist = eval('upl'+ str(x+1))
	pblist = eval('pb'+ str(x+1))
	groupline = 'group ' + plist + ', '
	for con in range(len(upllist)):
		# print(pblist[con])
		upl = upllist[con]
		cpline = 'copy ' + upl + ', ' + plist + '\n'
		groupline = groupline + upl + ' '
		pbout.write(pblist[con])
	outpml.write(groupline + "\n")
	pbout.close()
outpml.write("hide labels\n")


## Create poor constraints group in pml and write out value to summary files 
outcmx.write('open ' + outdir + outname + '_poor_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+3),'darkred'))
checkcons.write('### Low Support Constraints (final_poor_cons.pb) ###\n')
groupline = 'group poor cons ,'
for p in range(len(poorcons2)):
	groupline = groupline + poorcons[p] + ' '
	checkcons.write(poorcons2[p])
outpml.write(groupline + "\n")
checkcons.write('\n\n')
## Create long constraints group in pml and write out value to summary files 
outcmx.write('open ' + outdir +  outname + '_long_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+4),'aquamarine'))
checkcons.write('### Long Distance Constraints d >= 6.00 ###\n')
groupline = 'group long cons ,'
for l in range(len(longcons2)):
	checkcons.write(longcons2[l])
	groupline = groupline + longcons[l] + ' '
outpml.write(groupline + "\n")
checkcons.write('\n\n')

## Create short constraints group in pml and write out value to summary files 
outcmx.write('open ' + outdir + outname + '_short_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+5),'light coral'))
checkcons.write('### Short Distance Constraints d <= 3.00 ###\n')
groupline = 'group short cons ,'
for s in range(len(shortcons2)):
	checkcons.write(shortcons2[s])
	groupline = groupline + shortcons[s] + ' '
outpml.write(groupline + "\n")
checkcons.write('\n\n')
checkcons.close()
groupline = 'group viol peaks ,'
outcmx.write('open ' + outdir + outname + '_viol_peaks_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+6),'brown'))
for con in pmlpviols:
	groupline = groupline + con + ' '
outpml.write(groupline + "\n")

groupline = 'group viol upls ,'
outcmx.write('open ' + outdir + outname + '_viol_upl_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+7),'brown'))
for con in pmluviols:
	groupline = groupline + con + ' '
outcmx.write('open ' + outdir + 'hbond_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+8),'pink'))
outpml.write(groupline + "\n")
for uplfile in upls:
	x = x+8+1
	outcmx.write('open ' + outdir + uplfile.replace('.upl','_cons.pb') + '\n')
	outcmx.write('color #%s %s\n' %(str(x),'pink'))

selhbond = 'name hbond  #1.1:'
hbonsl = []
hbond = open(outdir + 'hbond_cons.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 0\n")
hbgroupline = 'group hbond , '
h = 1
for line in open('hbond.upl').readlines():
	cns = line.split()
	if "#" not in cns[0]:
		if cns[5] != 'H':
			h+=1 
			hbond.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], cns[2], cns[3],cns[5],'pink'))
			dist = 'distance hbond%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(h), pdbname, cns[0], cns[2], pdbname, cns[3], cns[5])
			tcolor = 'color pink, hbond' + str(h) + '\n'
			hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
			if cns[0] not in selhbond:
				selhbond = selhbond +'%s,' %(cns[0])
			if cns[3] not in selhbond:
				selhbond = selhbond +'%s,' %(cns[3])
		outpml.write(dist)
		outpml.write(tcolor)
outpml.write(hbgroupline + '\n')
selhbond = selhbond[:-1] + '@O,N\n'
outcmx.write(selhbond)
hbond.close()

outpml.write('color gray60, final\n')
outpml.write('split_states ' + pdbname + '\n')
for y in range(2,21,1):
	outcmx.write('match #1.%s to #1.1\n' %str(y))
	outpml.write('align %s_%04d, %s_0001\n' %(pdbname,y, pdbname))
outcmx.write('color :thr teal target a\n')
outcmx.write('color :thr teal target a\n')
outcmx.write('color :val orange  target a\n')
outcmx.write('color :leu indian red  target a\n')
outcmx.write('color :met purple  target a\n')
outcmx.write('color :ala forest green  target a\n')
outcmx.write('color :ile dodger blue  target a\n')
outcmx.write('color :phe slate blue target a\n')
outcmx.write('color :tyr orchid target a\n')
outcmx.write('color  byhetero target a\n')
outcmx.write('show :thr,met,ala,leu,val,ile,phe,tyr\n')
outcmx.write('hide H\n')
outcmx.write('show #1.1 @H\n')
outcmx.write('show #1.1@N target a\n')
outcmx.write('cartoon suppress false\n')
outcmx.write('label #1.1  text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\n')
outcmx.write('label ontop false\n')
outcmx.write('ui tool show "Side View"')
outpml.close()
outcmx.close()


