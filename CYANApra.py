### Mary Clay
import os
import sys
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
	name_pra.cxc
	name_pra.pml
	Pseudobond/Distance Groups from calculation:
		Each peak list 
		Poor constraints (SUP < 0.5)
		Short non intramolecular distance restraints d < 3.0
		Long distance restraints d > 6.0
		Peak Violations
		Manual restrain Violations (violated in 10 or more structures)
	Pseudobond/Distance Groups from manual restraints:
		input.upl
		hbond.upl 
	Anotated Constraints files 
''')
	exit()
## Dictionaries for assigning color to pml and cxc scripts for the peak list restraints
colors = ['white','raspberry','gold','forest','marine','purple','orange','cyan','pink','deepteal','gray']
colors2 = ['white','mediumvioletred','orange','forest','royalblue','purple','chocolate','cyan','pink','deepteal','gray']

cwd = os.getcwd() + '/'
outdir = cwd + 'post_cyana_ana/'
in_pdb = sys.argv[1]
fupl = sys.argv[2]
pdbname = in_pdb.split('.')[0]
fovw = fupl.replace('.upl','.ovw')
calc = cwd + 'CALC.cya'
outname = fupl.split('.')[0]

## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
	os.makedirs(outdir)

checkcons = open(outdir + outname + '_summary.txt','w')
checkcons.write('                         #peaks    upl Violations Assigned Ambiguous Unassigned\n')

## open the CALC.cya file to get the peaks list and additonal constraint files used in the calculation. 
cya_plists = [line.strip().replace('.peaks','-cycle7.peaks') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
print(manualongcons)
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
lols = [con for con in manualongcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualongcons if 'aco' in con]

print('                         #peaks    upl Violations Assigned Ambiguous Unassigned')
## Open Summary fiel and check the peak list files, upl, and ovw to determine the number of assignments and violations and write out the summary file 
## Creating the pseudobond files and group strings for rendering the constraints in chimera/pymol
for x in range(len(cya_plists)):
	plistn = cya_plists[x].replace('-cycle7.peaks','')
	exec("pb%s = open('%s','w')" %(str(x+1), outdir + outname + '_'+ plistn + '.pb'))
	pbout = eval('pb%s' %str(x+1))
	pbout.write("; halfbond = false\n; color = " + colors2[x+1] + "\n; radius = 0.1\n; dashes = 0\n")
	exec("group%s = '%s, '" %(str(x+1), 'group ' + cya_plists[x].replace('-cycle7.peaks','')))
	plist = cya_plists[x]
	upl = [line.strip() for line in open(fupl).readlines() if line.strip() and 'plist '+ str(x+1) in line]
	viol = [line.strip() for line in open(fovw).readlines() if line.strip() and 'list '+ str(x+1) in line]
	na,sa,aa = 0, 0, 0 
	cyplines = [line for line in open(os.getcwd() + '/' + plist).readlines() if line.strip() and line[0] != "#"]
	for line in cyplines:
		if 'e 0     0     0     0' in line: na+=1
		if len(line[0:8].strip()) !=0 and 'VC' in line: aa+=1
		if 'e 0     0     0     0' not in line and 'VC' not in line: sa+=1
	checkcons.write('%-25s%6d %6d %10d %8d %9d %10d\n' %(plist.replace('-cycle7.peaks',''),na+aa+sa,len(upl),len(viol),sa,aa,na))
	print('%-25s%6d %6d %10d %8d %9d %10d' %(plist.replace('-cycle7.peaks',''),na+aa+sa,len(upl),len(viol),sa,aa,na))
checkcons.write('\n\n')

outpml = open(outdir + fupl.replace('.upl','_pra.pml'),'w')
outpml.write('load '+ cwd + in_pdb+'\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('color gray60, all\n')
outcmx = open(outdir + fupl.replace('.upl','_pra.cxc'),'w')
outcmx.write('open '+ cwd + in_pdb+'\n')
outcmx.write('color #1 gray(150)\n')
pdbname = in_pdb.replace('.pdb','')

mns = len(cya_plists) + 8 + len(upls)
print(mns)

cmxphisel, cmxchisel, cmxphiviol, cmxchiviol = 'name phipsisel #%s:'%mns, 'name chisel #%s:'%(mns+1), 'name phipsiviol #%s:'%(mns+2), 'name chiviol #%s:'%(mns+3)
pmlphisel, pmlchisel, pmlphiviol, pmlchiviol = 'create phi-psi, %s_0001 and resi '%pdbname,'create chi, %s_0001 and resi ' %pdbname, 'create viol_phi-psi, %s_0001 and resi '%pdbname , 'create viol_chi, %s_0001 and resi '%pdbname


poorpbout = open(outdir + outname + '_poor_cons.pb','w')
poorpbout.write("; halfbond = false\n; color = darkred\n; radius = 0.2\n; dashes = 0\n")
longpbout = open(outdir + outname + '_long_cons.pb','w')
longpbout.write("; halfbond = false\n; color = aquamarine\n; radius = 0.2\n; dashes = 0\n")
shortpbout = open(outdir + outname + '_short_cons.pb','w')
shortpbout.write("; halfbond = false\n; color = light coral\n; radius = 0.2\n; dashes = 0\n")
pviolpbout = open(outdir + outname + '_viol_peaks_cons.pb','w')
pviolpbout.write("; halfbond = false\n; color = brown\n; radius = 0.2\n; dashes = 0\n")
uviolpbout = open(outdir + outname + '_viol_upls_cons.pb','w')
uviolpbout.write("; halfbond = false\n; color = brown\n; radius = 0.2\n; dashes = 0\n")

'''cns[0] = resi1 , cns[1] = resn1, cns[2]= atom1, cns[3] = resi2, cns[4]=resn2, cns[5] = atom2, cns[6] = dist'''
i = 1
Filtered = []
# name1 = line[10:15].strip()
# resn1 = line[16:20]
# resi1 = line[20:24].strip()
# name2 = line[27:32].strip()
# resn2 = line[33:37]
# resi2 = line[37:41].strip()
Upperdict, Lowerdict,= {}, {}
viol_peakscons,viol_uplscons= 'group viol_peaks, ', 'group viol_upls, '
finalupls = [["###Violated Restraints\n"],["###Poor/Low Support\n"],["###Long Distance Restraints (d >= 6.0)\n"],["###Short Distance Restraints (d <= 3.0)\n"],["###Good Restraints\n"]]
checkcons.write('### Violated Distance Constraints from %s \n' %(str(fovw)))
v = 0
phiv, chiv = [], []
for line in open(fovw).readlines():
	if line[4:9] == 'Upper' or line[4:9] == 'Lower':
		dviol = line.split()
		atom1 = dviol[1]
		if dviol[2]+dviol[1] in replacements.keys():
			atom1 = replacements[dviol[2]+dviol[1]]
		atom2 = dviol[5]
		if dviol[6]+dviol[5] in replacements.keys():
			atom2 = replacements[dviol[6]+dviol[5]]
		if dviol[9] >= '10':
			v+=1
			if 'peak' not in line:
				pbout = uviolpbout
				# grpout = 'viol_uplscons'
				grpstr = "uplviol"
				cons = '%4s %4s %-4s  %4s %4s %-4s  %6.2f\n' % (dviol[3],dviol[2],dviol[1],dviol[7],dviol[6],dviol[5],float(dviol[8]))
				cons2 = '%4s %4s %-4s  %4s %4s %-4s  %6.2f  # %s %s\n' % (dviol[3],dviol[2],dviol[1],dviol[7],dviol[6],dviol[5],float(dviol[8]), dviol[9], dviol[10])
				if line[4:9] == 'Upper':
					Upperdict[cons] = cons2
				if line[4:9] == 'lower':
					Lowerdict[cons] = cons2
			if 'peak' in line and 'QQ' not in line:
				pbout = pviolpbout
				# grpout = 'viol_peakscons'
				grpstr = "peakviol"
				for line2 in open(fupl).readlines():
					cns = line2.split()
					if cns[8] == line[90:].split()[1] and cns[10] == line[90:].split()[3] and cns[2] == dviol[1] and cns[5] == dviol[5]:
						checkcons.write(line2.replace('\n',' #Violated ' + line[50:88]+ '\n'))
						finalupls[0].append(line2)
						Filtered.append(line2)
			if ',' in atom1 and ',' not in atom2:
				pbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(dviol[3], atom1.split(',')[0], dviol[7],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(grpstr, str(v), pdbname, dviol[3], atom1.split(',')[0], pdbname, dviol[7], atom2))
				if 'peak' in line: viol_uplscons = viol_uplscons + "uplviol"+str(v) + ' '
				if 'peak' not in line: viol_peakscons = viol_peakscons + "peakviol"+str(v) + ' '
				v+=1
				pbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(dviol[3], atom1.split(',')[1], dviol[7],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(grpstr, str(v), pdbname, dviol[3], atom1.split(',')[1], pdbname, dviol[7], atom2))
				if 'peak' in line: viol_uplscons = viol_uplscons + "uplviol"+str(v) + ' '
				if 'peak' not in line: viol_peakscons = viol_peakscons + "peakviol"+str(v) + ' '
			if ',' in atom2 and ',' not in atom1:
				pbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(dviol[3], atom1, dviol[7],atom2.split(',')[0]))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(grpstr, str(v), pdbname, dviol[3], atom1, pdbname, dviol[7], atom2.split(',')[0]))
				if 'peak' in line: viol_uplscons = viol_uplscons + "uplviol"+str(v) + ' '
				if 'peak' not in line: viol_peakscons = viol_peakscons + "peakviol"+str(v) + ' '
				v+=1
				pbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(dviol[3], atom1, dviol[7],atom2.split(',')[1]))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(grpstr, str(v), pdbname, dviol[3], atom1, pdbname, dviol[7], atom2.split(',')[1]))
				if 'peak' in line: viol_uplscons = viol_uplscons + "uplviol"+str(v) + ' '
				if 'peak' not in line: viol_peakscons = viol_peakscons + "peakviol"+str(v) + ' '
			if ',' not in atom1 and ',' not in atom2:
				pbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(dviol[3], atom1, dviol[7],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(grpstr, str(v), pdbname, dviol[3], atom1, pdbname, dviol[7], atom2))
				if 'peak' in line: viol_uplscons = viol_uplscons + "uplviol"+str(v) + ' '
				if 'peak' not in line: viol_peakscons = viol_peakscons + "peakviol"+str(v) + ' '
	if line[4:9] == 'Angle':
		dang = line.split()
		if dang[1] == 'PHI' or dang[1] == 'PSI':
			if dang[3] not in phiv:
				phiv.append(dang[3])
				cmxphiviol = cmxphiviol + dang[3] + ','
				pmlphiviol = pmlphiviol + dang[3] + '+'
		if 'CHI' in dang[1]:
			if dang[3] not in chiv:
				chiv.append(dang[3])
				cmxchiviol = cmxchiviol + dang[3] + ','
				pmlchiviol = pmlchiviol + dang[3] + '+'


checkcons.write('\n\n')
finalupl,poorcons2, show, shortcons2,longcons2 = [],[],[],[],[]
poorcons, shortcons, longcons = 'group poor_cons, ', 'group short_cons, ', 'group long_cons, '

for line in open(fupl).readlines():
	if line not in Filtered:
		if 'QQ' in line.split()[2] or 'QQ' in line.split()[5]:
			pass 
		else:
			cns = line.split()
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
			if len(line.split()) < 12: ## do not have QU vaules, entries usually have SUP = 1.0 
				outpml.write('distance UPL%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
				pblist.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2))
				exec('group' + cns[10] + '=' + 'group' + cns[10] + '+ "UPL%s "' %str(i))
				Filtered.append(line)
				finalupls[4].append(line)
			if len(line.split()) >= 12:
				if float(cns[12]) > 0.5:
					if float(cns[6]) >= 6.0:
						longpbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2))
						longcons2.append(line)
						outpml.write('distance long%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
						longcons = longcons + 'long%s ' %str(i)
						Filtered.append(line)
						finalupls[2].append(line)
						# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
						# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
					if float(cns[6]) <= 3.00:
						if abs(int(cns[0])- int(cns[3])) > 1:
							# if atom1 != 'H' or atom2 != 'H':
								shortpbout.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'light coral'))
								shortcons2.append(line)
								outpml.write('distance short%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
								shortcons = shortcons + 'short%s ' %str(i)
								finalupls[3].append(line)
								Filtered.append(line)
								# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
								# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
					else:
						outpml.write('distance UPL%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
						pblist.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2))
						exec('group' + cns[10] + '=' + 'group' + cns[10] + '+ "UPL%s "' %str(i))
						Filtered.append(line)
						finalupls[4].append(line)
				if float(cns[12]) < 0.5:
					poorpbout.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2))
					poorcons2.append(line)
					outpml.write('distance poor%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
					poorcons = poorcons + 'poor%s ' %str(i)
					Filtered.append(line)
					finalupls[1].append(line)
					# show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
					# show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
poorpbout.close()
longpbout.close()
shortpbout.close()
pviolpbout.close()
uviolpbout.close()
#### Write out Poor/Low Support constraints to the summary file
checkcons.write('### Low Support Constraints (final_poor_cons.pb) ###\n')
for p in range(len(poorcons2)):
	checkcons.write(poorcons2[p])
checkcons.write('\n\n')
#### Write out Long Distance constraints to the summary file
checkcons.write('### Long Distance Constraints d >= 6.00 ###\n')
for l in range(len(longcons2)):
	checkcons.write(longcons2[l])
checkcons.write('\n\n')
#### Write out Short Distance constraints to the summary file
checkcons.write('### Short Distance Constraints d <= 3.00 ###\n')
for s in range(len(shortcons2)):
	checkcons.write(shortcons2[s])
checkcons.write('\n\n')
checkcons.close()

mn = 1
for x in range(len(cya_plists)):
	mn+=1
	pbout = eval('pb%s' %str(x+1))
	outcmx.write('open ' + outdir + outname + '_'+ cya_plists[x].replace('-cycle7.peaks','.pb\n'))
	outcmx.write('color #%s %s\n' %(str(mn),colors2[mn]))
	groupstr = eval('group' + str(x+1))
	outpml.write(groupstr + '\n')
	outpml.write('color %s, %s\n' %(colors[x+1],cya_plists[x].replace('-cycle7.peaks','')))

for (group, color) in [('poor','darkred'),('long','aquamarine'),('short', 'coral'),('viol_peaks', 'brown'),('viol_upls', 'brown')]:
	mn+=1
	outcmx.write('open ' + outdir + outname + '_' + group + '_cons.pb\n')
	outcmx.write('color #%s %s\n' %(str(mn),color))
	grpstr = eval(group + 'cons')
	outpml.write(grpstr + '\n')
	outpml.write('color %s, %s\n' %(color, group + 'cons'))
#### Write out the filtered upl list, which does not contain ambiguous (QQ) restraints 
#### and has sorted the restrints into 5 labeled catagories

filtered_upl = open(outdir + fupl.replace('.upl','4cns.upl'),'w')
for upllist in finalupls:
	for upl in upllist:
		filtered_upl.write(upl)
filtered_upl.close()

for uplfile in upls:
	newlines = []
	fin = open(uplfile,'r')
	for line in fin.readlines():
		if line in Upperdict.keys():
			newlines.append(Upperdict[line])
		if line not in Upperdict.keys():
			newlines.append(line)
	fout = open(outdir + uplfile,'w')
	fout.writelines(newlines)
fin.close()
fout.close()
for lolfile in lols:
	newlines = []
	fin = open(lolfile,'r')
	for line in fin.readlines():
		if line in Lowerdict.keys():
			newlines.append(Lowerdict[line])
		if line not in Lowerdict.keys():
			newlines.append(line)
	fout = open(outdir + lolfile,'w')
	fout.writelines(newlines)
fin.close()
fout.close()
u = 1
for uplfile in upls:
	fin = open(uplfile,'r')
	outpb = open(outdir + uplfile.replace('.upl','_cons.pb'),'w')
	pmlgroup = 'group %s, ' %(uplfile.replace('.upl',''))
	outpb.write("; halfbond = false\n; color = blue\n; radius = 0.1\n; dashes = 10\n")
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
				outpb.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1.split(',')[0], cns[3],atom2))
				outpb.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1.split(',')[1], cns[3],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1.split(',')[0], pdbname, cns[3], atom2))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1.split(',')[1], pdbname, cns[3], atom2))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
			if ',' in atom2 and ',' not in atom1:
				outpb.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2.split(',')[0]))
				outpb.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2.split(',')[1]))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2.split(',')[0]))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2.split(',')[2]))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
			if ',' not in atom1 and ',' not in atom2:
				outpb.write('#1.1:%s@%s #1.1:%s@%s\n' %(cns[0], atom1, cns[3],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
				pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
	outpml.write(pmlgroup + '\n')
	outpml.write('color blue,' + uplfile.replace('.upl','') + '\n')
	mn+=1
	outcmx.write('open ' + outdir + uplfile.replace('.upl','_cons.pb') + '\n')
	outcmx.write('color #%s %s\n' %(str(mn),'cyan'))

selhbond = 'name hbond  #1.1:'
hbonsl = []
hbond = open(outdir + 'hbond_cons.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 10\n")
hbgroupline = 'group hbond , '
h = 1
mn+=1
for line in open('hbond.upl').readlines():
	cns = line.split()
	if "#" not in cns[0]:
		if (cns[0],cns[3]) not in hbonsl:
			h+=1 
			hbonsl.append((cns[0],cns[3]))
			hbonsl.append((cns[3],cns[0]))
			hbond.write('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], cns[2], cns[3],cns[5],'pink'))
			outpml.write('distance hbond%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(h), pdbname, cns[0], cns[2].replace('H','N'), pdbname, cns[3], cns[5].replace('H','N')))
			hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
			if cns[0] not in selhbond:
				selhbond = selhbond +'%s,' %(cns[0])
			if cns[3] not in selhbond:
				selhbond = selhbond +'%s,' %(cns[3])
hbond.close()
outpml.write(hbgroupline + '\n')
outpml.write('color pink, hbond\n')
selhbond = selhbond[:-1] + '@O,N\n'
outcmx.write('open ' + outdir + 'hbond_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(mn),'pink'))
outcmx.write(selhbond)

phir, chir = [],[]
for angf in dihed:
	for line in open(angf).readlines():
		if '#' not in line and line.strip():
			ang = line.split()
			if ang[2] == 'PHI' or ang[2] == 'PSI':
				if ang[0] not in phir:
					phir.append(ang[0])
					cmxphisel = cmxphisel + ang[0] + ','
					pmlphisel = pmlphisel + ang[0] + '+'
			if 'CHI' in ang[2]:
				if ang[0] not in chir:
					chir.append(ang[0])
					cmxchisel = cmxchisel + ang[0] + ','
					pmlchisel = pmlchisel  + ang[0] + '+'


# outcmx.write('combine #1.1 modelId %s name angles\n' %mns)
outcmx.write('combine #1.1 modelId %s name phi-psi\n'%(mn+1))
outcmx.write(cmxphisel[:-1] + '\n')
outcmx.write('color phipsisel purple target ac\n')
outcmx.write('combine #1.1 modelId %s name chi\n'%(mn+2))
outcmx.write(cmxchisel[:-1] + '\n')
outcmx.write('color chisel cornflower blue target ac \n')
outcmx.write('combine #1.1 modelId %s name viol_phi-psi\n'%(mn+3))
outcmx.write(cmxphiviol[:-1] + '\n')
outcmx.write('color phipsiviol hot pink target ac \n')
outcmx.write('combine #1.1 modelId %s name viol_chi\n'%(mn+4))
outcmx.write(cmxchiviol[:-1] + '\n')
outcmx.write('color chiviol mediumvioletred target ac\n')

outpml.write(pmlphisel[:-1] + '\n')
outpml.write('color purple,phi-psi\n')
outpml.write(pmlchisel[:-1] + '\n')
outpml.write('color marin, chi\n')
outpml.write(pmlchiviol[:-1] + '\n')
outpml.write('color magenta, viol_phi-psi \n')
outpml.write(pmlphiviol[:-1] + '\n')
outpml.write('color red, viol_chi\n')
outpml.write("hide labels\n")
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


