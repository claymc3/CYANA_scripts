### Mary Clay
import os
import sys
import re
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
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH',
'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N'}

# ['ALAH','CYSH','ASPH','GLUH','PHEH','GLYH','HISH','ILEH','LYSH','LEUH','METH','ASNH','GLNH','ARGH','SERH','THRH','VALH','TRPH','TYRH']
if len(sys.argv)==1:
	print('''

Usage: 
	precya [pdb] [TALOS]

	precyana AF-FGFR3_KD.pdb ../TALOS2

Required Input:

	PDB			PDB to be used typically the final.pdb or pdb after CNS refinment
				If this is not located in current directory provide path
					CNS/refinePDB/r12_cya.pdb
	TALOS		Path to the TALOS results directory, scipt will find the predSS.tab
				file it needs 

OutPut:
	name_pra.cxc
	name_pra.pml
	colors cartoon representation based on talos secondary structure classification:
		Chemical Shift based Helical Navy
		Chemical Shift based Beta Strand Teal
		Chemical Shift based Loop Strand goldenrod
		Sequence based Helical royalblue
		Sequence based Beta Strand turquoise
		Sequence based Loop khaki
	Pseudobond/Distance Groups from upl files:
		model.upl
		hbond.upl
''')
	exit()
colors = ['orange','forest','royalblue','purple','chocolate','teal','gold','navy','darkturquoise','pink','cyan']

cwd = os.getcwd() + '/'
outdir = cwd + 'pre_cyana/'
in_pdb = sys.argv[1]
TALOSdir = sys.argv[-1]
print('in precyana')
print(TALOSdir)
talosSS = os.path.join(TALOSdir +'/predSS.tab')
pdbname = in_pdb.split('/')[-1].replace('.pdb','')
calc = cwd + 'CALC.cya'
outname = in_pdb.split('.')[0]

if '/' not in in_pdb:
	pdb_path = '../' + in_pdb
	pymol_pdb_path = './' + in_pdb

## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
	os.makedirs(outdir)

## open the CALC.cya file to get the peaks list and additonal constraint files used in the calculation. 
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
lols = [con for con in manualongcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualongcons if 'aco' in con]
mn = 1
outpml = open(outdir + 'CYANA_input.pml','w')
outpml.write('load '+pdb_path+'\n')
outpml.write('load '+pymol_pdb_path+'\n')
outpml.write('set_color teal = [0,127,127]\nset_color turquoise = [64,224,209]\nset_color goldenrod = [219,166,31]\nset_color palevioletred = [219,112,147]\nset_color orange = [255,165,0]\nset_color forest = [34,139,34]\nset_color royalblue = [65,105,225]\nset_color chocolate = [210,105,30]\nset_color purple = [128,0,128]\nset_color teal = [0,128,128]\nset_color gold = [255,215,0]\nset_color navy = [0,0,128]\nset_color darkturquoise = [0,206,209]\nset_color pink = [255,192,203]\nset_color cyan = [0,255,255]\nset_color paleturquoise = [175,238,238]\nset_color lightsalmon = [255,160,122]\nset_color khaki = [240,230,140]\nset_color yellowgreen = [154,205,50]\nset_color thistle = [216,191,216]\nset_color aquamarine = [127,255,212]\nset_color plum = [221,160,221]\nset_color lightpint = [255,182,193]\nset_color mediumvioletred = [199,21,133]\nset_color firebrick = [178,34,34]\nset_color lightcoral = [240,128,128]\nset_color deeppink = [255,20,147]\nset_color hotpink = [255,105,180]\nset_color purple = [128,0,128]\nset_color mediumpurple = [147,112,219]\nset_color navy = [0,0,128]\nset_color cornflowerblue = [100,149,237]\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('color gray60, all\n')
outcmx = open(outdir + 'CYANA_input.cxc','w')
outcmx.write('open {:} maxModels 1\nrename #{:} predSS\n'.format(pdb_path,mn))

mcount = 0
Hasprot = False
for line in open(in_pdb).readlines():
	if 'MODEL ' in line and 'REMARK' not in line:mcount+=1
	if line[0:4] == "ATOM" or line[0:4] == 'HETA':
		if line[12:16].strip() == 'H': Hasprot = True
if Hasprot == False:
	for aa in ['ALAH','CYSH','ASPH','GLUH','PHEH','GLYH','HISH','ILEH','LYSH','LEUH','METH','ASNH','GLNH','ARGH','SERH','THRH','VALH','TRPH','TYRH']:
		replacements[aa] = 'N'

if mcount > 2: 
	outpml.write('split_states ' + pdbname + '\n')
	pmln = '{:}_0001'.format(pdbname)
if mcount <= 1: 
	pmln = pdbname

### Color code secondar structure from TALOS analysis
outpml.write('create predSS, {:}\ncolor gray60,predSS\nhide sticks, predSS\n'.format(pmln))
outcmx.write('label #{:} text "{{0.label_one_letter_code}}{{0.number}}{{0.insertion_code}}"\nlabel ontop false\n'.format(mn))
CSHelix = 'name CSHelix #{:}:'.format(mn)
CSStrand = 'name CSStrand #{:}:'.format(mn)
CSLoop = 'name CSLoop #{:}:'.format(mn)
SeqHelix = 'name SeqHelix #{:}:'.format(mn)
SeqStrand = 'name SeqStrand #{:}:'.format(mn)
SeqLoop = 'name SeqLoop #{:}:'.format(mn)
talos_lines = [line.strip() for line in open(talosSS).readlines() if line.strip() and not re.search('[A-Z]', line[0])]
for line in talos_lines:
	if line.split()[-1] == 'H':CSHelix = CSHelix + line.split()[0] + ','
	if line.split()[-1] == 'E':CSStrand = CSStrand + line.split()[0] + ','
	if line.split()[-1] == 'L':CSLoop = CSLoop + line.split()[0] + ','
	if line.split()[-1] == 'h':SeqHelix = SeqHelix + line.split()[0] + ','
	if line.split()[-1] == 'e':SeqStrand = SeqStrand + line.split()[0] + ','
	if line.split()[-1] == 'c':SeqLoop = SeqLoop + line.split()[0] + ','

if CSHelix[-1] != ":":
	outcmx.write(CSHelix[:-1]+ '\ncolor CSHelix navy target c\n')
	outpml.write('color navy, predSS and resi ' + CSHelix[CSHelix.index(':')+1:-1].replace(',','+') + '\n')
if CSStrand[-1]!= ":":
	outcmx.write(CSStrand[:-1]+ '\ncolor CSStrand teal target c\n')
	outpml.write('color teal, predSS and resi ' + CSStrand[CSStrand.index(':')+1:-1].replace(',','+') + '\n')
if CSLoop[-1] != ":":
	outcmx.write(CSLoop[:-1]+ '\ncolor CSLoop goldenrod target c\n')
	outpml.write('color goldenrod, predSS and resi ' + CSLoop[CSLoop.index(':')+1:-1].replace(',','+') + '\n')
if SeqHelix[-1] != ":":
	outcmx.write(SeqHelix[:-1]+ '\ncolor SeqHelix royal blue target c\n')
	outpml.write('color royalblue, predSS and resi ' + SeqHelix[SeqHelix.index(':')+1:-1].replace(',','+') + '\n')
if SeqStrand[-1] != ":":
	outcmx.write(SeqStrand[:-1]+ '\ncolor SeqStrand turquoise target c\n')
	outpml.write('color turquoise, predSS and resi ' + SeqStrand[SeqStrand.index(':')+1:-1].replace(',','+') + '\n')
if SeqLoop[-1] != ":":
	outcmx.write(SeqLoop[:-1]+ '\ncolor SeqLoop khaki target c\n')
	outpml.write('color khaki, predSS and resi ' + SeqLoop[SeqLoop.index(':')+1:-1].replace(',','+') + '\n')
### Read in hbond and render the pseudo bonds on the predSS model in chimera 
selhbond = 'name hbond  #{:}:'.format(mn)
sidehbond = 'name shbond #{}:'.format(mn)
hbonsl = []
hbond = open(outdir + 'hbond.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 5\n")
hbgroupline = 'group hbond , '
h = 0
for line in open('hbond.upl').readlines():
	cns = line.split()
	if line.split():
		if "#" not in cns[0]:
			if (cns[0],cns[3]) not in hbonsl:
				h+=1 
				hbonsl.append((cns[0],cns[3]))
				#hbonsl.append((cns[3],cns[0]))
				hbond.write('#{:}:{:}@{:} #{:}:{:}@{:}\n'.format(mn, cns[0], cns[2], mn, cns[3],cns[5]))
				outpml.write('distance hbond{:}, {:} and resi {:} and name {:}, {:} and resi {:} and name {:}\n'.format(str(h), pmln, cns[0], cns[2].replace('H','N'), pmln, cns[3], cns[5].replace('H','N')))
				hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
				if cns[0] not in selhbond:
					selhbond = selhbond +'{:},'.format(cns[0])
				if cns[3] not in selhbond:
					selhbond = selhbond +'{:},'.format(cns[3])
				if cns[2] not in ['O','H','N']: sidehbond = sidehbond +'{:},'.format(cns[0])
				if cns[5] not in ['O','H','N']: sidehbond = sidehbond +'{:},'.format(cns[3])
hbond.close()
outpml.write(hbgroupline + '\n')
outcmx.write('open '+ pdb_path +' maxModels 1\n')
cmxn = '#2'
outcmx.write('color {:} gray(150)\n'.format(cmxn))
mn+=1
### Read in the dihed.aco file and color residues that have defined phi/psi angles purple, and defined chi angles cornflower blue
cmxphisel, cmxchisel = 'name phipsisel #{:}:'.format(mn), 'name chisel #{:}:'.format(mn)
pmlphisel, pmlchisel = 'color purple, phi-psi and resi ','color navy, chi and resi '
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
outcmx.write('open {:} maxModels 1\nrename #{:} dihed\nhide #{:} target a\n color #{:} gray(150)\n'.format(pdb_path,mn,mn,mn))
outcmx.write(cmxphisel[:-1] + '\n')
outcmx.write('color phipsisel purple target ac\n')
outpml.write('create phi-psi, {:}\ncolor gray60, phi-psi\nhide sticks, phi-psi\n'.format(pmln))
outpml.write(pmlphisel[:-1] + '\n')
sidelist = []
u = 0
x = -1
mn+=1
for uplfile in upls:
	x+=1
	fin = open(uplfile,'r')
	outpb = open(outdir + uplfile.replace('.upl','.pb'),'w')
	outpb.write("; halfbond = false\n; color = {:}\n; radius = 0.15\n; dashes = 10\n".format(colors[x]))
	gid = 'group {:}, '.format(uplfile.replace('.upl',''))
	for line in fin.readlines():
		if line.split():
			cns = line.split()
			if "#" not in cns[0]:
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
				atoms2 = atom2.split(',')
				atoms1 = atom1.split(',')
				for atom1 in atoms1:
					for atom2 in atoms2:
						u+=1
						outpb.write('{:}:{:}@{:} {:}:{:}@{:}\n'.format(cmxn, cns[0], atom1, cmxn, cns[3],atom2))
						outpml.write('distance {:}{:}, {:} and resi {:} and name {:}, {:} and resi {:} and name {:}\n'.format(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1, pmln, cns[3], atom2))
						gid = gid + "{:}{:} ".format(uplfile.replace('.upl',''),str(u))
				if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
					sidelist.append(cns[0])
				if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
					sidelist.append(cns[3])
	mn+=1
	if re.search(uplfile.replace('.upl','')+'[0-9]*',gid):
		outpml.write(gid + '\n')
		outpml.write('color {:}, {:}\n'.format(colors[x],uplfile.replace('.upl','')))
		outcmx.write('open ' + outdir + uplfile.replace('.upl','.pb') + '\n')
		# print('{:} {:}'.format(mn,uplfile.replace('.upl','.pb')))
		outcmx.write('color #{:} {:}\n'.format(str(mn),colors[x]))
	if not re.search(uplfile.replace('.upl','')+'[0-9]*',gid):
		os.remove(outdir + uplfile.replace('.upl','.pb'))

sidechains = 'show #2:'
for res in sidelist:
	sidechains = sidechains + res + ','
outcmx.write(sidechains[:-1] + '\n')
outpml.write('show sticks, {:} and resi {:}\n'.format(pdbname, sidechains[sidechains.index(':')+1:-1].replace(',','+')))

outpml.write('color pink, hbond\n')
selhbond = selhbond[:-1] + '@O,N\n'
sidehbond = sidehbond[:-1] + '\n'
outcmx.write('open ' + outdir + 'hbond.pb\n')
mn+=1
outcmx.write('color #{:} {:}\n'.format(str(mn),'pink'))
outcmx.write(selhbond)
outcmx.write(sidehbond)

if len(cmxchisel[:-1]) > 18:
	outcmx.write(cmxchisel[:-1] + '\n''color chisel navy target a \n')
	outcmx.write('show chisel\n')
	outpml.write('create chi, {:}\ncolor gray60, chi\nshow sticks, chi\n'.format(pmln))
	outpml.write(pmlchisel[:-1] + '\n')
outpml.write('color gray60, {:}\n'.format(pdbname))
outpml.write('show sticks, {:} and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\n hide sticks, elem H\nhide sticks, name N+C\n'.format(pdbname))
outpml.write('color paleturquoise, {:} and resn ILE\ncolor lightsalmon, {:} and resn LEU\ncolor khaki, {:} and resn VAL\ncolor yellowgreen, {:} and resn ALA\ncolor thistle, {:} and resn MET\ncolor aquamarine, {:} and resn THR\ncolor lightpink, {:} and resn TYR\ncolor plum, {:} and resn PHE\n'.format(pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname))
outpml.write('color gold, elem S\ncolor red, elem O\ncolor blue, elem N\n')
outcmx.write('cartoon suppress false\nlabel {:}  text "{{0.label_one_letter_code}}{{0.number}}{{0.insertion_code}}"\nlabel ontop false\n'.format(cmxn))
outcmx.write('color #2:ile paleturquoise target a\ncolor #2:leu lightsalmon  target a\ncolor #2:val khaki target a\ncolor #2:ala yellowgreen  target a\ncolor #2:met thistle target a\ncolor #2:thr aquamarine target a\ncolor #2:phe plum target a\ncolor #2:tyr lightpink target a\n')
outcmx.write('color  byhetero target a\n')
outcmx.write('show #2:thr,met,ala,leu,val,ile,phe,tyr\n')
outcmx.write('hide H\nshow {:}@C,O,N,H target a\n'.format(cmxn))
outcmx.write('hide #{:}@C,O,N,H\n'.format(mn+1))
outcmx.write('show shbond  target a\nhide H\nshow hbond target a\n')
outcmx.write('ui tool show "Side View"\n')
outpml.write("hide labels\n")
for y in range(2,21,1):
	outpml.write('delete {:}_{:04d}\n'.format(pdbname,y))
outpml.close()
outcmx.close()


