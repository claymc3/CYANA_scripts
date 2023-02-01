### Mary Clay
import os
import sys
import re
import pandas as pd
import numpy as np
import GetDihe as Dihed
import noa_analysis as noaa

## HIST is his protonated on NE2
## HIS+ is his protonated on both ND1 and NE2 
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB',
'CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB',
'ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB',
'GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG',
'PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2',
'GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA',
'HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1',
'HISTHA':'CA','HISTHB2':'CB','HISTHB3':'CB','HISTQB':'CB','HISTHD1':'ND1','HISTHE2':'NE2','HISTHD2':'CD2','HISTHE1':'CE1',
'HIS+HA':'CA','HIS+HB2':'CB','HIS+HB3':'CB','HIS+QB':'CB','HIS+HD1':'ND1','HISE+E2':'NE2','HIS+HD2':'CD2','HIS+HE1':'CE1',
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
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH',
'PTRHA':'CA','PTRHB2':'CB','PTRHB3':'CB','PTRQB':'CB','PTRQD':'CD1,CD2','PTRQE':'CE1,CE2','PTRHD1':'CD1','PTRHE1':'CE1','PTRHE2':'CE2','PTRHD2':'CD2','PTRHH':'OH'}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P","cPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S'}

Ambiguous = {'ARGQB':'HB2,HB3','ARGQG':'HG2,HG3','ARGQG':'HG2,HG3','ARGQH1':'HH11,HH12','ARGQH2':'HH21,HH22','ASNQB':'HB2,HB3','ASNQD2':'HD21,HD22','ASPQB':'HB2,HB3','CYSQB':'HB2,HB3','CYSSQB':'HB2,HB3','GLNQB':'HB2,HB3','GLNQG':'HG2,HG3','GLNQE2':'HE21,HE22','GLUQB':'HB2,HB3','GLUQG':'HG2,HG3','GLYQA':'HA2,HA3','HISQB':'HB2,HB3','HISTQB':'HB2,HB3','ILEQG1':'HG12,HG13','LEUQB':'HB2,HB3','LEUQQD':'QD1,QD2','LYSQB':'HB2,HB3','LYSQG':'HG2,HG3','LYSQD':'HD2,HD3','LYSQE':'HE2,HE3','METQB':'HB2,HB3','METQG':'HG2,HG3','PHEQB':'HB2,HB3','PHEQG':'HG2,HG3','PHEQD':'HD2,HD3','PHEQE':'HE1,HE2','PROQB':'HB2,HB3','PROQG':'HG2,HG3','PROQD':'HD2,HD3','SERQB':'HB2,HB3','TRPQB':'HB2,HB3','TYRQB':'HB2,HB3','TYRQG':'HG2,HG3','TYRQD':'HD2,HD3','TYRQE':'HE1,HE2','VALQQG':'QG1,QG2'}

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
colors = ['white','palevioletred','orange','forest','royalblue','purple','chocolate','teal','gold','navy','darkturquoise','pink','cyan']


cwd = os.getcwd() + '/'
outdir = cwd + 'post_cyana_analysis/'
noadir = cwd  +'noa_analysis/'
in_pdb = sys.argv[1]
fupl = sys.argv[2]
pdbname = in_pdb.split('/')[-1].split('.')[0]
fovw = fupl.replace('.upl','.ovw')
calc = cwd + 'CALC.cya'
outname = fupl.split('.')[0]
init = cwd + 'init.cya'
# print(open(init).readlines()[0].strip().split(':=')[-1])
seq = [line.strip().split() for line in open(cwd + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0]]
Seqdict = {}
Sequence, ASequence = [],[]
for resn,resi in seq:
	Seqdict[resi] = AAA_dict[resn] + resi
	ASequence.append(AAA_dict[resn] + resi)
	Sequence.append(resi)
upldf = pd.DataFrame(index = Sequence, columns=['cya','long','viol','input','viol input','vdihed'])
upldf['cya'] = np.zeros(len(Sequence))
upldf['long'] = np.zeros(len(Sequence))
upldf['viol'] = np.zeros(len(Sequence))
upldf['input'] = np.zeros(len(Sequence))
upldf['viol input'] = np.zeros(len(Sequence))
## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
	os.makedirs(outdir)
if not os.path.exists(outdir +'pseudobonds/'):
	os.makedirs(outdir +'pseudobonds/')
if not os.path.exists(noadir):
	os.makedirs(noadir)
checkcons = open(outdir + outname + '_summary.txt','w')
## open the CALC.cya file to get the peaks list and additonal constraint files used in the calculation. 
cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
lengths = []
for plist in cya_plists:
	lengths.append(len(plist))
pad = ''
for x in range(max(lengths)):
	pad = pad + ' '
checkcons.write('{:}                                         Assignments \n{:}#peaks   upl  Viol Unique  Multiple  Unused  None  Diagonal Increased upl\n'.format(pad,pad))
manualcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualcons if 'upl' in con and 'hbond' not in con]
hbonds = [con for con in manualcons if 'upl' in con and 'hbond' in con]
lols = [con for con in manualcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualcons if 'aco' in con if con != 'inital.aco']
noa = cwd + 'cycle7.noa'
noalines = open(noa).readlines()
print('{:}                                         Assignments \n{:}#peaks   upl  Viol Unique  Multiple  Unused  None  Diagonal Increased upl'.format(pad,pad))
## Open Summary file and check the peak list files, upl, and ovw to determine the number of assignments and violations and write out the summary file 
## Creating the pseudobond files and group strings for rendering the constraints in chimera/pymol
tpeak,tsingle,tamb,tnotused,tnota,tdia,tincr,tupl,tviol  = 0, 0, 0, 0, 0, 0, 0, 0, 0
for x in range(len(cya_plists)):
	plistn = cya_plists[x]
	exec("pb{:} = open('{:}','w')".format(str(x+1), outdir +'pseudobonds/' + plistn + '.pb'))
	pbout = eval('pb{:}'.format(str(x+1)))
	pbout.write("; halfbond = false\n; color = " + colors[x+1] + "\n; radius = 0.1\n; dashes = 0\n")
	exec("group{:} = '{:}, '".format(str(x+1), 'group ' + cya_plists[x]))
	plist = cya_plists[x]
	upl = [line.strip() for line in open(fupl).readlines() if line.strip() and 'plist '+ str(x+1) in line]
	tupl = tupl + len(upl)
	viol = [line.strip() for line in open(fovw).readlines() if line.strip() and 'list '+ str(x+1) in line]
	tviol = tviol + len(viol)
	peak,single,amb,notused,nota,dia,incr = 0, 0, 0, 0, 0, 0, 0
	for y in range(len(noalines)):
		line = noalines[y]
		if ' ' + plistn in line and 'out of' in noalines[y+1]:
			peak+=1
			tpeak+=1
			if '0 out of 0' in noalines[y+1]:
				nota+=1
				tnota+=1
			if '0 out of' in noalines[y+1] and '0 out of 0' not in noalines[y+1]:
				notused+=1
				tnotused+=1
			if '1 out of' in noalines[y+1] and 'diagonal' not in line:
				single+= 1
				tsingle+=1
			if noalines[y+1].strip().split()[0] > '1' and 'diagonal' not in line:
				amb+=1
				tamb+=1
			if 'increased' in line:
				incr+=1 
				tincr+=1
			if 'diagonal' in line and '0 out of' not in noalines[y+1]:
				dia+=1
				tdia+=1
	linepad = pad[len(plist):]
	checkcons.write("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%\n".format(plist,linepad,peak,len(upl),len(viol), single,amb,notused,nota,dia,incr, 100*((single+amb+dia)/peak)))
	print("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%".format(plist,linepad, peak,len(upl),len(viol),single,amb,notused,nota,dia,incr, 100*((single+amb+dia)/peak)))
checkcons.write("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%\n".format('Total',pad[5:],tpeak,tupl,tviol,tsingle,tamb,tnotused,tnota,tdia,tincr,100*((tsingle+ tamb+ tdia)/tpeak)))
print("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%".format('Total',pad[5:],tpeak,tupl,tviol,tsingle,tamb,tnotused,tnota,tdia,tincr, 100*((tsingle+ tamb+ tdia)/tpeak)))
checkcons.write('\n\n')

outpml = open(outdir + fupl.replace('.upl','_pra.pml'),'w')
outpml.write('load '+ cwd + in_pdb+'\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('set_color palevioletred = [219,112,147]\nset_color orange = [255,165,0]\nset_color forest = [34,139,34]\nset_color royalblue = [65,105,225]\nset_color chocolate = [210,105,30]\nset_color purple = [128,0,128]\nset_color teal = [0,128,128]\nset_color gold = [255,215,0]\nset_color navy = [0,0,128]\nset_color darkturquoise = [0,206,209]\nset_color pink = [255,192,203]\nset_color cyan = [0,255,255]\nset_color paleturquoise = [175,238,238]\nset_color lightsalmon = [255,160,122]\nset_color khaki = [240,230,140]\nset_color yellowgreen = [154,205,50]\nset_color thistle = [216,191,216]\nset_color aquamarine = [127,255,212]\nset_color plum = [221,160,221]\nset_color lightpint = [255,182,193]\nset_color mediumvioletred = [199,21,133]\nset_color firebrick = [178,34,34]\nset_color lightcoral = [240,128,128]\nset_color deeppink = [255,20,147]\nset_color hotpink = [255,105,180]\nset_color purple = [128,0,128]\nset_color mediumpurple = [147,112,219]\nset_color navy = [0,0,128]\nset_color cornflowerblue = [100,149,237]\n')
outpml.write('color gray60, all\n')
outpml.write('show sticks, {:} and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\n hide sticks, elem H\nhide sticks, name N+C\n'.format(pdbname))
outpml.write('color paleturquoise, {:} and resn ILE\ncolor lightsalmon, {:} and resn LEU\ncolor khaki, {:} and resn VAL\ncolor yellowgreen, {:} and resn ALA\ncolor thistle, {:} and resn MET\ncolor aquamarine, {:} and resn THR\ncolor lightpink, {:} and resn TYR\ncolor plum, {:} and resn PHE\n'.format(pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname))
outpml.write('color gold, elem S\ncolor red, elem O\ncolor blue, elem N\n')
outpml.write('split_states ' + pdbname + '\n')
outcmx = open(outdir + fupl.replace('.upl','_pra.cxc'),'w')
outcmx.write('open '+ cwd + in_pdb+'\n')
outcmx.write('color #1 gray(150)\n')
outcmx.write('match #1.2-20 to #1.1\n')
outcmx.write('color #1:ile paleturquoise target a\ncolor #1:leu lightsalmon  target a\ncolor #1:val khaki target a\ncolor #1:ala yellowgreen  target a\ncolor #1:met thistle target a\ncolor #1:thr aquamarine target a\ncolor #1:phe plum target a\ncolor #1:tyr lightpink target a\n')
outcmx.write('show #1:thr,met,ala,leu,val,ile,phe,tyr\nname meyfside #1:thr,met,ala,leu,val,ile,phe,tyr\ncartoon suppress false\n')
outcmx.write('label #1.1 text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\n''label ontop false\n')
outcmx.write('ui tool show "Side View"\n#ui mousemode right distance\n')

angmn = len(cya_plists) + 8 + len(upls)
cmxphisel, cmxchisel, cmxphiviol, cmxchiviol = 'name phipsisel #{:}:'.format(angmn), 'name chisel #{:}:'.format(angmn), 'name phipsiviol #{:}:'.format(angmn), 'name chiviol #{:}:'.format(angmn)
pmlphisel, pmlchisel, pmlphiviol, pmlchiviol = 'color purple, phi-psi and resi ','color navy, chi and resi ', 'color mediumpurple, viol_phi-psi and resi ', 'color cornflowerblue, viol_chi and resi '
poorpbout = open(outdir +'pseudobonds/' + outname + '_poor_cons.pb','w')
poorpbout.write("; halfbond = false\n; color = mediumvioletred\n; radius = 0.1\n; dashes = 0\n")
longpbout = open(outdir +'pseudobonds/' + outname + '_long_cons.pb','w')
longpbout.write("; halfbond = false\n; color = firebrick\n; radius = 0.1\n; dashes = 0\n")
shortpbout = open(outdir +'pseudobonds/' + outname + '_short_cons.pb','w')
shortpbout.write("; halfbond = false\n; color = lightcoral\n; radius = 0.1\n; dashes = 0\n")
pviolpbout = open(outdir +'pseudobonds/' + outname + '_viol_peaks_cons.pb','w')
pviolpbout.write("; halfbond = false\n; color = deeppink\n; radius = 0.1\n; dashes = 0\n")
uviolpbout = open(outdir +'pseudobonds/' + outname + '_viol_upls_cons.pb','w')
uviolpbout.write("; halfbond = false\n; color = hotpink\n; radius = 0.1\n; dashes = 0\n")

### Go through the final overview file and extract information about violated distance and angle restraints 

i = 1
Filtered = []
violdict, Upperdict, Lowerdict = {}, {}, {}
viol_peakscons,viol_uplscons= 'group viol_peaks, ', 'group viol_upls, '
finalupls = [["###Violated Restraints\n"],["###Poor/Low Support\n"],["###Long Distance Restraints (d >= 6.0)\n"],["###Short Distance Restraints (d <= 3.0)\n"],["###Good Restraints\n"]]
checkcons.write('### Violated Distance Constraints from {:} \n'.format(str(fovw)))
v = 0
phiv, chiv, violpeaks = [], [], []
for line in open(fovw).readlines():
	if line[4:9] == 'Upper' or line[4:9] == 'Lower':
		dviol = line.split()
		atom1 = dviol[1]
		if dviol[2]+dviol[1] in replacements.keys():
			atom1 = replacements[dviol[2]+dviol[1]]
		atom2 = dviol[5]
		if dviol[6]+dviol[5] in replacements.keys():
			atom2 = replacements[dviol[6]+dviol[5]]
		atoms1 = atom1.split(',')
		atoms2 = atom2.split(',')
		cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[2]],dviol[3],dviol[1],AAA_dict[dviol[6]],dviol[7],dviol[5])
		cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[6]],dviol[7],dviol[5],AAA_dict[dviol[2]],dviol[3],dviol[1])
		if int(dviol[9]) >= 10:
			v+=1
			if 'peak' not in line:
				pbout = uviolpbout
				grpstr = "uplviol"
				if line[4:9] == 'Upper':
					outline = ' #viol in {:} by +{:}\n'.format(dviol[9],dviol[10])
					upldf.loc[dviol[7],'viol input'] = upldf.loc[dviol[7],'viol input'] + 1
					upldf.loc[dviol[3],'viol input'] = upldf.loc[dviol[3],'viol input'] + 1
					Upperdict[cons1] = outline
					Upperdict[cons2] = outline
				if line[4:9] == 'lower':
					outline = ' #viol in {:} by -{:}\n'.format(dviol[9],dviol[10])
					Lowerdict[cons1] = outline
					Lowerdict[cons2] = outline
			if 'peak' in line and 'list' in line:
				pbout = pviolpbout
				grpstr = "peakviol"
				outline = ' #viol in {:} by {:}\n'.format(dviol[9],dviol[10])
				violdict[cons1] = outline
				violdict[cons2] = outline
				for line2 in open(fupl).readlines():
					cns = line2.split()
					if cns[8] == line[90:].split()[1] and cns[10] == line[90:].split()[3] and cns[2] == dviol[1] and cns[5] == dviol[5]:
						violpeaks.append(line2.replace('\n',' #Violated ' + line[44:88]+ '\n'))
						if cns[0] != cns[3] and '#SUP' in line2:
							upldf.loc[cns[0],'viol'] = upldf.loc[cns[0],'viol'] + 1
							upldf.loc[cns[3],'viol'] = upldf.loc[cns[3],'viol'] + 1
							finalupls[0].append(line2)
							Filtered.append(line2)
			for atom1 in atoms1:
				for atom2 in atoms2: 
					v+=1
					pbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(dviol[3], atom1, dviol[7],atom2))
					outpml.write('distance viol{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(v), pdbname, dviol[3], atom1, pdbname, dviol[7], atom2))
					if 'peak' in line: viol_uplscons = viol_uplscons + "viol"+str(v) + ' '
					if 'peak' not in line: viol_peakscons = viol_peakscons + "viol"+str(v) + ' '
	if line[4:9] == 'Angle':
		dang = line.split()
		angle = upldf.loc[dang[3],'vdihed']
		if pd.isna(angle): angle = ''
		angle = angle + '{:} {:}\n'.format(dang[1],dang[7])
		upldf.loc[dang[3],'vdihed'] = angle
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
violpeaks = sorted(violpeaks, key = lambda x: (x.split()[10],x.split()[8]))
for viol in violpeaks:
	checkcons.write(viol)
checkcons.write('\n\n')
usedupls,qupldict, upldict, upldict2 = {}, {}, {}, []
finalupl,poorcons2, show, shortcons2,longcons2,sidelist = [],[],[],[],[],[]
poorcons, shortcons, longcons = 'group poor, ', 'group short, ', 'group long, '
for line in open(fupl).readlines():
	if line.split() and '#SUP' in line:
		cns = line.split()
		atom1 = cns[2]
		atom2 = cns[5]
		cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],atom1,AAA_dict[cns[4]],cns[3],atom2)
		cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],atom2,AAA_dict[cns[1]],cns[0],atom1)
		if cons1 not in upldict.keys():
			upldict[cons1] = "{:3.2f}".format(float(cns[6]))
		if cons2 not in upldict.keys():
			upldict[cons2] = "{:3.2f}".format(float(cns[6]))
		upldict2.append("{:} peak {:} from {:}".format(cons1,cns[8],cya_plists[int(cns[10])-1]))
		upldict2.append("{:} peak {:} from {:}".format(cons2,cns[8],cya_plists[int(cns[10])-1]))

for line in open(fupl).readlines():
	if line not in Filtered:
		if '#SUP' not in line: ## exclude ambiguous restraints
			pass 
		else:
			cns = line.split()
			if cns[0] == cns[3]: ## exclude intramolecular restraints
				pass
			if cns[0] != cns[3]:
				pblist = eval('pb' + cns[10])
				atom1 = cns[2]
				atom2 = cns[5]
				if cns[1]+cns[2] in replacements.keys():
					atom1 = atom1.replace(cns[2], replacements[cns[1]+cns[2]])
				if cns[4]+cns[5] in replacements.keys():
					atom2 = atom2.replace(cns[5], replacements[cns[4]+cns[5]])
				atoms1 = atom1.split(',')
				atoms2 = atom2.split(',')
				upldf.loc[cns[0],'cya'] = upldf.loc[cns[0],'cya'] + 1
				upldf.loc[cns[3],'cya'] = upldf.loc[cns[3],'cya'] + 1
				cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
				cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5],AAA_dict[cns[1]],cns[0],cns[2])
				outline = ' #{:3.2f} #peak {:} #plist {:}\n'.format(float(cns[6]),cns[8],cns[10])
				usedupls[cons1] = outline
				usedupls[cons2] = outline
				for atom1 in atoms1:
					if atom1 == 'H': atom1 = 'N'
					for atom2 in atoms2:
						if atom2 == 'H': atom2 = 'N'
						cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],atom1,AAA_dict[cns[4]],cns[3],atom2)
						cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],atom2,AAA_dict[cns[1]],cns[0],atom1)
						usedupls[cons1] = outline
						usedupls[cons2] = outline
				if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
					sidelist.append(cns[0])
				if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
					sidelist.append(cns[3])
				if float(cns[12]) < 0.5:
					for atom1 in atoms1:
						for atom2 in atoms2:
							i+=1
							poorpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
							poorcons2.append(line)
							outpml.write('distance poor{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
							poorcons = poorcons + 'poor{:} '.format(i)
					Filtered.append(line)
					finalupls[1].append(line)
				if float(cns[12]) > 0.5:
					if float(cns[6]) >= 6.0:
						upldf.loc[cns[0],'long'] = upldf.loc[cns[0],'long'] + 1
						upldf.loc[cns[3],'long'] = upldf.loc[cns[3],'long'] + 1
						cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
						cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5],AAA_dict[cns[1]],cns[0],cns[2])
						qupldict[cons1] = " #long distance\n"
						qupldict[cons2] = " #long distance\n"
						for atom1 in atoms1:
							for atom2 in atoms2:
								i+=1
								longpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
								longcons2.append(line)
								outpml.write('distance long{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
								longcons = longcons + 'long{:} '.format(i)
						Filtered.append(line)
						finalupls[2].append(line)
					if float(cns[6]) <= 3.00:
						if abs(int(cns[0])- int(cns[3])) > 1:
							for atom1 in atoms1:
								for atom2 in atoms2:
									i+=1
									shortpbout.write('#1.1:{:}@{:} #1.1:{:}@{:} {:}\n'.format(cns[0], atom1, cns[3],atom2, 'light coral'))
									shortcons2.append(line)
									outpml.write('distance short{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
									shortcons = shortcons + 'short{:} '.format(i)
							finalupls[3].append(line)
							Filtered.append(line)
					if float(cns[6]) > 3.00 and float(cns[6]) < 6.0:
						for atom1 in atoms1:
							for atom2 in atoms2:
								i+=1
								outpml.write('distance UPL{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
								pblist.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
								exec('group' + cns[10] + '=' + 'group' + cns[10] + '+ "UPL{:} "'.format(i))
						Filtered.append(line)
						finalupls[4].append(line)
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
	pbout = eval('pb{:}'.format(str(x+1)))
	outcmx.write('open ' + outdir +'pseudobonds/' + cya_plists[x] +'.pb\n')
	outcmx.write('color #{:} {:}\n'.format(str(mn),colors[mn]))
	groupstr = eval('group' + str(x+1))
	outpml.write(groupstr + '\n')
	outpml.write('color {:}, {:}\n'.format(colors[x+1],cya_plists[x]))

for (group, color) in [('poor','mediumvioletred'),('long','firebrick'),('short', 'lightcoral'),('viol_peaks', 'deeppink'),('viol_upls', 'hotpink')]:
	mn+=1
	outcmx.write('open ' + outdir +'pseudobonds/' + outname + '_' + group + '_cons.pb\n')
	outcmx.write('color #{:} {:}\n'.format(str(mn),color))
	grpstr = eval(group + 'cons')
	outpml.write(grpstr + '\n')
	outpml.write('color {:}, {:}\n'.format(color, group))
#### Write out the filtered upl list, which does not contain ambiguous (QQ) restraints 
#### and has sorted the restrints into 5 labeled catagories

filtered_upl = open(outdir + fupl.replace('.upl','4cns.upl'),'w')
for upllist in finalupls:
	for upl in upllist:
		filtered_upl.write(upl)
filtered_upl.close()

u = 1
for uplfile in upls:
	fin = open(uplfile,'r')
	outpb = open(outdir +'pseudobonds/' + uplfile.replace('.upl','.pb'),'w')
	pmlgroup = 'group {:}, '.format(uplfile.replace('.upl',''))
	outpb.write("; halfbond = false\n; color = blue\n; radius = 0.15\n; dashes = 10\n")
	for line in fin.readlines():
		cns = line.split()
		if line.strip() and "#" not in cns[0]:
			atom1 = cns[2]
			atom2 = cns[5]
			if cns[1]+cns[2] in replacements.keys():
				atom1 = atom1.replace(cns[2], replacements[cns[1]+cns[2]])
			if cns[4]+cns[5] in replacements.keys():
				atom2 = atom2.replace(cns[5], replacements[cns[4]+cns[5]])
			atoms2 = atom2.split(',')
			atoms1 = atom1.split(',')
			upldf.loc[cns[0],'input'] = upldf.loc[cns[0],'input'] + 1
			upldf.loc[cns[3],'input'] = upldf.loc[cns[3],'input'] + 1
			for atom1 in atoms1:
				for atom2 in atoms2:
					u+=1
					outpb.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
					outpml.write('distance {:}{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
					pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
			if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
				sidelist.append(cns[0])
			if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
				sidelist.append(cns[3])
	outpml.write(pmlgroup + '\n')
	outpml.write('color cyan,' + uplfile.replace('.upl','') + '\n')
	mn+=1
	outcmx.write('open ' + outdir +'pseudobonds/' + uplfile.replace('.upl','.pb') + '\n')
	outcmx.write('color #{:} {:}\n'.format(str(mn),'cyan'))
sidechains = 'show #1:'
for res in sidelist:
	sidechains = sidechains + res + ','
outcmx.write(sidechains[:-1] + '\n')
outcmx.write('hide H\n''show #1.1@H,N target a\n')

selhbond = 'name hbond  #1.1:'
hbonsl = []
hbond = open(outdir +'pseudobonds/' + 'hbond.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 10\n")
hbgroupline = 'group hbond , '
h = 1
mn+=1
for hbondf in hbonds:
	for line in open(hbondf).readlines():
		cns = line.split()
		if line.strip() and "#" not in cns[0]:
			if (cns[0],cns[3]) not in hbonsl:
				h+=1 
				hbonsl.append((cns[0],cns[3]))
				# hbonsl.append((cns[3],cns[0]))
				hbond.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], cns[2], cns[3],cns[5]))
				outpml.write('distance hbond{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(h), pdbname, cns[0], cns[2].replace('H','N'), pdbname, cns[3], cns[5].replace('H','N')))
				hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
				if cns[0] not in selhbond:
					selhbond = selhbond +'{:},'.format(cns[0])
				if cns[3] not in selhbond:
					selhbond = selhbond +'{:},'.format(cns[3])
hbond.close()
outpml.write(hbgroupline + '\n')
outpml.write('color pink, hbond\n')
selhbond = selhbond[:-1] + '@O,N\nshow hbond target a\n'
outcmx.write('open ' + outdir +'pseudobonds/' + 'hbond.pb\n')
outcmx.write('color #{:} {:}\n'.format(str(mn),'pink'))
outcmx.write(selhbond)
upls.extend(hbonds)
## Examin input upl files and update lines for entries which are violated 10 or more times. 
for uplfile in upls:
	newlines = []
	violupl = open(outdir + uplfile.replace('.upl','_viol.upl'),'w')
	matchedupl = open(outdir + uplfile.replace('.upl','_found.upl'),'w')
	for line in open(uplfile).readlines():
		newline = ''
		if line.strip():
			if '#' not in line.split()[0]:
				cns = line.split()
				cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
				if cons in Upperdict.keys():
					newline = line.replace('\n', Upperdict[cons])
					violupl.write(newline)
				if cons in usedupls.keys():
					if len(newline) > 1: newline = newline.replace('\n', usedupls[cons])
					if len(newline) < 1: newline = line.replace('\n', usedupls[cons])
					matchedline = line.replace('\n', usedupls[cons])
					matchedupl.write(matchedline)
		if len(newline) < 1:
			newline = line
		newlines.append(newline)
	fout = open(outdir + uplfile,'w')
	fout.writelines(newlines)
	fout.close()
	violupl.close()
	matchedupl.close()
for lolfile in lols:
	newlines = []
	for line in open(lolfile).readlines():
		newline = ''
		if line.strip():
			if '#' not in line.split()[0]:
				cns = line.split()
				cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
				if cons in Lowerdict.keys():
					newline = line.replace('\n', Lowerdict[cons])
		if len(newline) < 1:
			newline = line
		newlines.append(newline)
	fout = open(outdir + lolfile,'w')
	fout.writelines(newlines)
	fout.close()
print('finished finding upls')
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
for y in range(2,21,1):
	outpml.write('align {:}_{:04d}, {:}_0001\n'.format(pdbname,y, pdbname))
outcmx.write('open '+ cwd + in_pdb+ ' maxModels 1\nrename #{:} angles\nhide #{:} target a\n color #{:} gray(150)\n'.format(angmn,angmn,angmn))
outcmx.write(cmxphisel[:-1] + '\n')
outcmx.write('color phipsisel purple target c\n')
if cmxchisel[-1] != ':':
	outcmx.write(cmxchisel[:-1] + '\n')
	outcmx.write('color chisel navy target a \n')
	outcmx.write('show chisel target a\n')
if cmxphiviol[-1] != ':':
	outcmx.write(cmxphiviol[:-1] + '\n')
	outcmx.write('color phipsiviol mediumpurple target c \n')
if cmxchiviol[-1] != ':':
	outcmx.write(cmxchiviol[:-1] + '\n')
	outcmx.write('color chiviol cornflower blue target a\n')
	outcmx.write('show chiviol target a\n')
outcmx.write('hide #{:}@H*,N,O target a\ncolor  byhetero target a\n'.format(angmn))
outpml.write('hide everything, {:}\n'.format(pdbname))
outpml.write('create phi-psi, {:}_0001\ncolor gray60,phi-psi\nhide sticks, phi-psi\n'.format(pdbname))
outpml.write(pmlphisel[:-1] + '\n')
outpml.write('create chi, {:}_0001\ncolor gray60, chi\nhide sticks, chi\n'.format(pdbname))
outpml.write(pmlchisel[:-1] + '\n')
outpml.write('create viol_phi-psi, {:}_0001\ncolor gray60, viol_phi-psi\nhide sticks, viol_phi-psi\n'.format(pdbname))
outpml.write(pmlchiviol[:-1] + '\n')
outpml.write('create viol_chi, {:}_0001\ncolor gray60, viol_chi\nhide sticks, viol_chi\n'.format(pdbname))
outpml.write(pmlphiviol[:-1] + '\n')
outpml.write("hide labels\n")
outpml.close()
outcmx.close()

# ---------------------------------------------------------------------------
# Run the cycle7.noa analysis and generate peak list identifying peaks as 
# unused, no assingment, and questionable

noaa.analize_noa(cwd, noadir, calc, noa, Seqdict, violdict, qupldict,upldict, pad, upldict2)
# ---------------------------------------------------------------------------
# Run the GetDihed.py to determine phi, psi, chi1 and chi2 and plot them
# for all 20 structures

Dihed.extract(in_pdb, ASequence, outdir)
print('finished plotting dihedrals')

# ---------------------------------------------------------------------------
# Plot the upl counts for each residue 
# 

import matplotlib as mpl 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplcursors
import numpy as np
from matplotlib.widgets import Slider
print(pdbname)
pdf=PdfPages(outdir + '{:}_upl_overview.pdf'.format(pdbname))
nsubplots = round(len(Sequence)/25,0)
if round(len(Sequence)/25,1) - nsubplots > 0.0:
	nsubplots = nsubplots + 1
if nsubplots == 0: nsubplots = 1
fig_height = 3.0 * nsubplots
if fig_height <= 2.0: 
	fig_height = 3.0
fig=plt.figure(figsize=(10.0,fig_height))
spi = 0
for i in range(0,len(upldf.index.tolist()),25):
	spi = spi + 1
	temp, resid = [], []
	z = i
	for y in range(25):
		temp.append(upldf.index.tolist()[z])
		resid.append(ASequence[z])
		z = z + 1 
		if z == len(upldf.index.tolist()):break
	index = np.arange(i,z,1)
	dfp = upldf.reindex(temp)
	width = 0.18
	ax = fig.add_subplot(int(nsubplots),1,spi)
	ax.bar(index, dfp['cya'], width, color = '#9acd32', ecolor='none', label='CYANA UPL')
	ax.bar(index + width, dfp['long'], width, color = '#800080', edgecolor='none', label='long UPL')
	ax.bar(index + 2 * width, dfp['viol'], width, color = '#ffa500', edgecolor='none', label='Violated UPL')
	ax.bar(index + 3 * width, dfp['input'], width, color = '#6495ed', edgecolor='none', label='Input UPL')
	ax.bar(index + 4 * width, dfp['viol input'], width, color = '#db7093', edgecolor='none', label='Violate Input UPL')
	angelidx, angelh = [], []
	for n in range(len(index)):
		res = temp[n]
		angle = dfp.loc[res,'vdihed']
		if not pd.isna(angle):
			angelidx.append(index[n]+2*width)
			angelh.append(max(upldf['cya']))
			text = ax.text(index[n]+0.1, max(upldf['cya']), angle,ha='center', va='top',fontsize=6)
	ax.bar(angelidx, angelh, 0.9, color='gray',alpha = 0.5, zorder = 0.0,edgecolor='none' )
	ax.set_xticks(index +2*width, labels=resid)
	ax.set_ylim([0,max(upldf['cya'])])
	ax.set_xlim([(min(index) - 1.0 * width), (max(index) + 5 * width)])
	ax.tick_params(axis='x', labelrotation = 90)
	ax.legend(loc='best', frameon=False, markerscale=0.000001)
	ax.set_ylabel('Number of UPL Entries')
	ax.set_xlabel('Residue')
	# cursors = mplcursors.cursor(hover=True)
	plt.tight_layout(pad = 0.4, w_pad = 0.4, h_pad = 0.4)
pdf.savefig()
# plt.show()
plt.close()
pdf.close()
print('finished')

