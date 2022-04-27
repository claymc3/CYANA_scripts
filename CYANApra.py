import os
import sys
import glob
import pandas as pd
replacements ={'ALAH':'N','ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB',
'CYSH':'N','CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB',
'ASPH':'N','ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB',
'GLUH':'N','GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG',
'PHEH':'N','PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2',
'GLYH':'N','GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA',
'HISH':'N','HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1',
'ILEH':'N','ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEQD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1',
'LYSH':'N','LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ',
'LEUH':'N','LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUQQD':'CD2,CD1',
'METH':'N','METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHG2':'CG','METHG3':'CG','METQG':'CG','METQE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE',
'ASNH':'N','ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD2':'ND2',
'PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD',
'GLNH':'N','GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNQE2':'NE2',
'ARGH':'N','ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2',
'SERH':'N','SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG',
'THRH':'N','THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2',
'VALH':'N','VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALQQG':'CG1,CG2',
'TRPH':'N','TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2',
'TYRH':'N','TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH'}




calc = os.getcwd() + '/CALC.cya'
fovw = os.getcwd() + '/final.ovw'

in_pdb = sys.argv[1]
fupl = sys.argv[2]
pdbname = in_pdb.split('.')[0]

summary = pd.DataFrame(columns=['#peaks', 'upl', 'Violations', 'Assigned', 'Ambiguous', 'Unassigned' ])
cya_plists = [line.strip().replace('.peaks','-cycle7.peaks') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')

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
	summary.loc[plist.replace('-cycle7.peaks',''),'#peaks'] = na+aa+sa
	summary.loc[plist.replace('-cycle7.peaks',''),'upl'] = len(upl)
	summary.loc[plist.replace('-cycle7.peaks',''),'Violations'] = len(viol)
	summary.loc[plist.replace('-cycle7.peaks',''),'Unassigned'] = na
	summary.loc[plist.replace('-cycle7.peaks',''),'Ambiguous'] = aa
	summary.loc[plist.replace('-cycle7.peaks',''),'Assigned'] = sa
print(summary)

colors = ['white','raspberry','gold','forest','marine','purple','orange','cyan','pink','deepteal','gray']
colors2 = ['white','mediumvioletred','orange','forest','royalblue','purple','chocolate','cyan','pink','deepteal','gray']

outpml = open(fupl.replace('.upl','_dist.pml'),'w')
outpml.write('load '+ in_pdb+'\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('color gray60, all\n')
outcmx = open(fupl.replace('.upl','_dist.cxc'),'w')
outcmx.write('open '+ in_pdb+'\n')
pdbname = in_pdb.replace('.pdb','')
'''cns[0] = resi1 , cns[1] = resn1, cns[2]= atom1, cns[3] = resi2, cns[4]=resn2, cns[5] = atom2, cns[6] = dist'''
i = 1
finalupl = []
finalupl,badcons,badcons2, show, scons, scons2, lcons,lcons2 = [],[],[],[],[],[],[],[]
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
		upllist.append('UPL' + str(i))
		dist = 'distance UPL %s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2)
		tcolor = 'color ' + colors[int(cns[10])]  + ', UPL' + str(i) + '\n'
		pblist.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, colors2[int(cns[10])]))
		i+=1
		if len(cns) >= 12:
			if float(cns[12]) < 0.5:
				badcons.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'darkred'))
				badcons2.append(line)
				show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
				show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
			if float(cns[12]) > 0.5:
				finalupl.append(line)
		if float(cns[6]) >= 6.0:
			lcons.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'cyan'))
			lcons2.append(line)
			show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
			show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
		if float(cns[6]) <= 3.00:
			if cns[0] != cns[3]:
				if atom1 != 'H' or atom2 != 'H':
					scons.append('#1.1:%s@%s #1.1:%s@%s %s\n' %(cns[0], atom1, cns[3],atom2, 'cyan'))
					scons2.append(line)
					show.append('show #1.1:%s@%s target a\n' %(cns[0], atom1))
					show.append('show #1.1:%s@%s target a\n' %(cns[3], atom2))
					finalupl.append(line)
					#finalupl.append(line.replace(cns[6], '3.60').replace('\n','  ## adjusted to VDW dist\n'))
		outpml.write(dist)
		outpml.write(tcolor)

filtered_upl = open(fupl.replace('.upl','4cns.upl'),'w')
for upl in finalupl:
	filtered_upl.write(upl)
filtered_upl.close()
# print(badcons2)
for x in range(len(cya_plists)):
	plist = cya_plists[x].replace('-cycle7.peaks','')
	outname = fupl.split('.')[0]
	outcmx.write('open ' + outname + '_'+ plist + '.pb\n')
	outcmx.write('color #%s %s\n' %(str(x+2),colors2[x+1]))
	pbout = open(outname + '_'+ plist + '.pb','w')
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

checkcons = open(outname + '_check_cons.txt','w')
outcmx.write('open ' + outname + '_poor_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+3),'darkred'))
pbout = open(outname + '_poor_cons.pb','w')
pbout.write("; halfbond = false\n")
pbout.write("; color = darkred\n")
pbout.write("; radius = 0.2\n")
pbout.write("; dashes = 0\n")
checkcons.write('### Low Support Constraints (final_poor_cons.pb) ###\n')
for con in badcons:
	pbout.write(con)
pbout.close()
for con in badcons2:
	checkcons.write(con)

outcmx.write('open ' + outname + '_long_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+4),'cyan'))
pbout = open(outname + '_long_cons.pb','w')
pbout.write("; halfbond = false\n")
pbout.write("; color = darkred\n")
pbout.write("; radius = 0.2\n")
pbout.write("; dashes = 0\n")
for con in lcons:
	pbout.write(con)
pbout.close()
checkcons.write('### Long Distance Constraints d >= 6.00 ###\n')
for con in lcons2:
	checkcons.write(con)

outcmx.write('open ' + outname + '_short_cons.pb\n')
outcmx.write('color #%s %s\n' %(str(x+5),'orange'))
pbout = open(outname + '_short_cons.pb','w')
pbout.write("; halfbond = false\n")
pbout.write("; color = darkred\n")
pbout.write("; radius = 0.2\n")
pbout.write("; dashes = 0\n")
for con in scons:
	pbout.write(con)
pbout.close()
checkcons.write('### Short Distance Constraints d <= 3.00 ###\n')
for con in scons2:
	checkcons.write(con)
checkcons.close()

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
outcmx.write('cartoon suppress false\n')
outcmx.write('label #1.1  text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\n')
outcmx.write('label ontop false\n')

for con in show:
	outcmx.write(con)
outpml.close()
outcmx.close()


fovw = os.getcwd() + '/final.ovw'
upls = glob.glob(os.getcwd() + '/*.upl')
upls = [upl for upl in upls if 'cycle' not in upl]
lols = glob.glob(os.getcwd() + '/*.lol')
# name1 = line[10:15].strip()
# resn1 = line[16:19]
# resi1 = line[20:24].strip()
# name2 = line[27:32].strip()
# resn2 = line[33:36]
# resi2 = line[37:41].strip()
Upperdict, Lowerdict,= {}, {}
for line in open('final.ovw').readlines():
	if line[4:9] == 'Upper' or line[4:9] == 'Lower':
		if 'peak' not in line:
			cons = '%4s %s  %-3s   %4s %s  %-3s   %6.2f\n' % (line[20:24].strip(),line[16:19],line[10:15].strip(),line[37:41].strip(),line[33:36],line[27:32].strip(),float(line[44:48]))
			cons2 = '## %4s %s  %-3s   %4s %s  %-3s   %6.2f  #%s %s\n' % (line[20:24].strip(),line[16:19],line[10:15].strip(),line[37:41].strip(),line[33:36],line[27:32].strip(),float(line[44:48]), line[50:52], line[62:66])
			if line[4:9] == 'Upper':
				Upperdict[cons] = cons2
			if line[4:9] == 'lower':
				Lowerdict[cons] = cons2

# print(Upperdict.keys())
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
