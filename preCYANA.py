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
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH'}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',

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
	TALOS		Path to the TALOS results directory, scipt will finde the predSS.tab
				file it needs 
				Chemical Shift based Helical Navy
				Chemical Shift based Betta Strand Teal
				Chemical Shift based Loop Strand goldenrod
				Sequence based Helical royalblue
				Sequence based Betta Strand turquoise
				Sequence based Loop khaki

OutPut:
	name_pra.cxc
	name_pra.pml
	Pseudobond/Distance Groups from manual restraints:
		input.upl
		hbond.upl 
''')
	exit()
colors = ['mediumvioletred','orange','forest','royalblue','purple','chocolate','cyan','palevioletred','deepteal','gold','navy','darkcyan']

cwd = os.getcwd() + '/'
outdir = cwd + 'pre_cyana/'
in_pdb = sys.argv[1]
TALOSdir = sys.argv[2]
talosSS = os.path.join(TALOSdir +'/predSS.tab')
pdbname = in_pdb.split('.')[0]
calc = cwd + 'CALC.cya'
outname = in_pdb.split('.')[0]

## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
	os.makedirs(outdir)



## open the CALC.cya file to get the peaks list and additonal constraint files used in the calculation. 
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
lols = [con for con in manualongcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualongcons if 'aco' in con]

outpml = open(outdir + 'CYANA_input.pml','w')
outpml.write('load '+ cwd + in_pdb+'\n')
outpml.write('set_color navy = [0,0,128]\nset_color royalblue = [65,105,225]\nset_color darkcyan = [0,139,139]\nset_color turquoise = [64,224,208]\nset_color goldenrod = [218,165,32]\nset_color khaki = [240,230,140]\nset_color mediumvioletred = [199,21,133]\nset_color gold = [255,215,0]\nset_color cornflowerblue = [100,149,237]\nset_color teal = [0,128,128]\nset_color palevioletred = [219,112,147]\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('color gray60, all\n')
outcmx = open(outdir + 'CYANA_input.cxc','w')
outcmx.write('open '+ cwd + in_pdb+'\n')
outcmx.write('color #1 gray(150)\n')
pdbname = in_pdb.replace('.pdb','')

mn = 1
mcount = 0
Hasprot = False
for line in open(in_pdb).readlines():
	if 'MODEL ' in line:mcount+=1
	if line[0:4] == "ATOM" or line[0:4] == 'HETA':
		if line[12:16].strip() == 'H': Hasprot = True
if Hasprot == False:
	for aa in ['ALAH','CYSH','ASPH','GLUH','PHEH','GLYH','HISH','ILEH','LYSH','LEUH','METH','ASNH','GLNH','ARGH','SERH','THRH','VALH','TRPH','TYRH']:
		replacements[aa] = 'N'

if mcount > 2: 
	cmxn = '#1.1'
	pmln = '%s_0001' %pdbname
if mcount <= 1: 
	cmxn = '#1'
	pmln = pdbname


u = 0
x = -1
for uplfile in upls:
	mn+=1
	x+=1
	fin = open(uplfile,'r')
	NNpb = open(outdir + uplfile.replace('.upl','_NN.pb'),'w')
	NNpb.write("; halfbond = false\n; color = %s\n; radius = 0.1\n; dashes = 0\n" %colors[x])
	NCpb = open(outdir + uplfile.replace('.upl','_NC.pb'),'w')
	NCpb.write("; halfbond = false\n; color = %s\n; radius = 0.1\n; dashes = 0\n" %colors[x+1])
	CCpb = open(outdir + uplfile.replace('.upl','_CC.pb'),'w')
	CCpb.write("; halfbond = false\n; color = %s\n; radius = 0.1\n; dashes = 0\n" %colors[x+2])
	groupNN = 'group %s, ' %(uplfile.replace('.upl','_NN'))
	groupNC = 'group %s, ' %(uplfile.replace('.upl','_NC'))
	groupCC = 'group %s, ' %(uplfile.replace('.upl','_CC'))
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
			if atom1[0] == 'N' and atom2[0] == 'N':
				outpb = NNpb
				gid = 'NN'
			if atom1[0] == 'N' and atom2[0] == 'C':
				outpb = NCpb
				gid = 'NC'
			if atom1[0] == 'C' and atom2[0] == 'C':
				outpb = CCpb
				gid = 'CC'
			if ',' in atom1 and ',' not in atom2:
				outpb.write('%s:%s@%s %s:%s@%s\n' %(cmxn, cns[0], atom1.split(',')[0], cmxn,  cns[3],atom2))
				outpb.write('%s:%s@%s %s:%s@%s\n' %(cmxn, cns[0], atom1.split(',')[1], cmxn, cns[3],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1.split(',')[0], pmln, cns[3], atom2))
				exec('group' + gid + '=' + 'group' + gid + '+"%s%s "' %(uplfile.replace('.upl',''),str(u)))
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1.split(',')[1], pmln, cns[3], atom2))
				exec('group' + gid + '=' + 'group' + gid + '+"%s%s "' %(uplfile.replace('.upl',''),str(u)))
			if ',' in atom2 and ',' not in atom1:
				outpb.write('%s:%s@%s %s:%s@%s\n' %(cmxn, cns[0], atom1, cmxn, cns[3],atom2.split(',')[0]))
				outpb.write('%s:%s@%s %s:%s@%s\n' %(cmxn, cns[0], atom1, cmxn, cns[3],atom2.split(',')[1]))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1, pmln, cns[3], atom2.split(',')[0]))
				exec('group' + gid + '=' + 'group' + gid + '+"%s%s "' %(uplfile.replace('.upl',''),str(u)))
				u+=1
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1, pmln, cns[3], atom2.split(',')[2]))
				outpml.write('color pink, %s%s\n'%(uplfile.replace('.upl',''),str(u)))
				exec('group' + gid + '=' + 'group' + gid + '+"%s%s "' %(uplfile.replace('.upl',''),str(u)))
			if ',' not in atom1 and ',' not in atom2:
				outpb.write('%s:%s@%s %s:%s@%s\n' %(cmxn, cns[0], atom1, cmxn, cns[3],atom2))
				outpml.write('distance %s%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(uplfile.replace('.upl',''),str(u), pmln, cns[0], atom1, pmln, cns[3], atom2))
				exec('group' + gid + '=' + 'group' + gid + '+"%s%s "' %(uplfile.replace('.upl',''),str(u)))

	if re.search(uplfile.replace('.upl','')+'[0-9]*',groupNN):
		outpml.write(groupNN + '\n')
		outpml.write('color %s, %s\n' %(colors[x],uplfile.replace('.upl','_NN')))
		outcmx.write('open ' + outdir + uplfile.replace('.upl','_NN.pb') + '\n')
		print('%s %s' %(mn,uplfile.replace('.upl','_NN.pb')))
		outcmx.write('color #%s %s\n' %(str(mn),colors[x]))
	if not re.search(uplfile.replace('.upl','')+'[0-9]*',groupNN):
		os.remove(outdir + uplfile.replace('.upl','_NN.pb'))
	if re.search(uplfile.replace('.upl','')+'[0-9]*',groupNC):
		mn+=1
		x+=1
		outpml.write(groupNC + '\n')
		outpml.write('color %s, %s\n' %(colors[x+1],uplfile.replace('.upl','_NC')))
		outcmx.write('open ' + outdir + uplfile.replace('.upl','_NC.pb') + '\n')
		print('%s %s' %(mn,uplfile.replace('.upl','_NC.pb')))
		outcmx.write('color #%s %s\n' %(str(mn),colors[x+1]))
	if not re.search(uplfile.replace('.upl','')+'[0-9]*',groupNC):
		os.remove(outdir + uplfile.replace('.upl','_NC.pb'))
	if re.search(uplfile.replace('.upl','')+'[0-9]*',groupCC):
		mn+=1
		x+=1
		outpml.write(groupCC + '\n')
		outpml.write('color %s, %s\n' %(colors[x+2],uplfile.replace('.upl','_CC')))
		outcmx.write('open ' + outdir + uplfile.replace('.upl','_CC.pb') + '\n')
		print('%s %s' %(mn,uplfile.replace('.upl','_CC.pb')))
		outcmx.write('color #%s %s\n' %(str(mn),colors[x+2]))
	if not re.search(uplfile.replace('.upl','')+'[0-9]*',groupCC):
		os.remove(outdir + uplfile.replace('.upl','_CC.pb'))



### Color code secondar structure from TALOS analysis
outpml.write('create predSS, %s\ncolor gray60,predSS\nhide sticks, predSS\n' %pmln)
outcmx.write('open '+ cwd + in_pdb+'\n')
mn+=1
if mcount > 2: 
	hbcmxn = '#%s.1'%mn
	outcmx.write('rename #%s predSS\n' %(mn))
	outcmx.write('hide #%s.2-20 target ac\n'%mn)
if mcount <= 1: 
	hbcmxn = '#%s' %mn
	outcmx.write('rename #%s predSS\n' %(mn))
outcmx.write('label %s  text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\nlabel ontop false\n' %hbcmxn)
CSHelix = 'name CSHelix %s:' %hbcmxn
CSStrand = 'name CSStrand %s:' %hbcmxn
CSLoop = 'name CSLoop %s:' %hbcmxn
SeqHelix = 'name SeqHelix %s:' %hbcmxn
SeqStrand = 'name SeqStrand %s:' %hbcmxn
SeqLoop = 'name SeqLoop %s:' %hbcmxn
talos_lines = [line.strip() for line in open(talosSS).readlines() if line.strip() and re.search(' *([0-9]*) [A-Z]', line)]
for line in talos_lines:
	if line.split()[-1] == 'H':CSHelix = CSHelix + line.split()[0] + ','
	if line.split()[-1] == 'E':CSStrand = CSStrand + line.split()[0] + ','
	if line.split()[-1] == 'L':CSLoop = CSLoop + line.split()[0] + ','
	if line.split()[-1] == 'h':SeqHelix = SeqHelix + line.split()[0] + ','
	if line.split()[-1] == 'e':SeqStrand = SeqStrand + line.split()[0] + ','
	if line.split()[-1] == 'c':SeqLoop = SeqLoop + line.split()[0] + ','

outcmx.write(CSHelix[:-1]+ '\ncolor CSHelix navy target c\n')
outcmx.write(SeqHelix[:-1]+ '\ncolor SeqHelix royal blue target c\n')
outcmx.write(CSStrand[:-1]+ '\ncolor CSStrand teal target c\n')
outcmx.write(SeqStrand[:-1]+ '\ncolor SeqStrand turquoise target c\n')
outcmx.write(CSLoop[:-1]+ '\ncolor CSLoop goldenrod target c\n')
outcmx.write(SeqLoop[:-1]+ '\ncolor SeqLoop khaki target c\n')

outpml.write('color navy, predSS and resi ' + CSHelix[CSHelix.index(':'):-1].replace(',','+') + '\n')
outpml.write('color royalblue, predSS and resi ' + SeqHelix[SeqHelix.index(':'):-1].replace(',','+') + '\n')
outpml.write('color teal, predSS and resi ' + CSStrand[CSStrand.index(':'):-1].replace(',','+') + '\n')
outpml.write('color turquoise, predSS and resi ' + SeqStrand[SeqStrand.index(':'):-1].replace(',','+') + '\n')
outpml.write('color goldenrod, predSS and resi ' + CSLoop[CSLoop.index(':'):-1].replace(',','+') + '\n')
outpml.write('color khaki, predSS and resi ' + SeqLoop[SeqLoop.index(':'):-1].replace(',','+') + '\n')

### Read in hbond and render the pseudo bonds on the predSS model in chimera 

selhbond = 'name hbond  #%s:'%mn
hbonsl = []
hbond = open(outdir + 'hbond_cons.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 5\n")
hbgroupline = 'group hbond , '
h = 0
for line in open('hbond.upl').readlines():
	cns = line.split()
	if "#" not in cns[0]:
		if (cns[0],cns[3]) not in hbonsl:
			h+=1 
			hbonsl.append((cns[0],cns[3]))
			hbonsl.append((cns[3],cns[0]))
			hbond.write('%s:%s@%s %s:%s@%s\n' %(hbcmxn, cns[0], cns[2], hbcmxn, cns[3],cns[5]))
			outpml.write('distance hbond%s, %s and resi %s and name %s, %s and resi %s and name %s\n' %(str(h), pmln, cns[0], cns[2].replace('H','N'), pmln, cns[3], cns[5].replace('H','N')))
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
mn+=1
outcmx.write('color #%s %s\n' %(str(mn),'pink'))
outcmx.write(selhbond)


### Read in the dihed.aco file and color residues that have defined phi/psi angles purple, and defined chi angles cornflower blue
cmxphisel, cmxchisel = 'name phipsisel #%s:'%(mn+1), 'name chisel #%s:'%(mn+2)
pmlphisel, pmlchisel = 'color purple, phi-psi and resi ','color cornflowerblue, chi and resi '
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

outcmx.write('combine %s modelId %s name phi-psi\n'%(cmxn, mn+1))
outcmx.write(cmxphisel[:-1] + '\n')
outcmx.write('color phipsisel purple target ac\n')
outpml.write('create phi-psi, %s\ncolor gray60, phi-psi\nhide sticks, phi-psi\n' %pmln)

if len(cmxchisel[:-1]) > 18:
	outcmx.write('combine %s modelId %s name chi\n'%(cmxn, mn+2))
	outcmx.write(cmxchisel[:-1] + '\n')
	outcmx.write('color chisel cornflower blue target ac \n')
	outpml.write('create chi, %s\ncolor gray60, chi\nhide sticks, chi\n' %pmln)
	outpml.write(pmlphisel[:-1] + '\n')



outpml.write('color gray60, %s\n' %pdbname)
outpml.write('show sticks, %s and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\n hide sticks, elem H\nhide sticks, name N+C\n' %pdbname)
outpml.write('color violetpurple, %s and resn MET\ncolor smudge, %s and resn ALA\ncolor marine, %s and resn ILE\ncolor magenta, %s and resn LEU\ncolor orange, %s and resn VAL\ncolor gold, %s and resn THR\ncolor pink, %s and resn TYR\ncolor slate, %s and resn PHE\n' %(pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname))
outpml.write('color gold, elem S\ncolor red, elem O\ncolor blue, elem N\n')

outcmx.write('color #1:thr teal target a\ncolor #1:val orange  target a\ncolor #1:leu indian red  target a\ncolor #1:met purple  target a\ncolor #1:ala forest green  target a\ncolor #1:ile dodger blue  target a\ncolor #1:phe slate blue target a\ncolor #1:tyr orchid target a\n')
outcmx.write('color  byhetero target a\n')
outcmx.write('show #1:thr,met,ala,leu,val,ile,phe,tyr\n')
outcmx.write('hide H\nshow %s@N,H target a\n' %cmxn)
outcmx.write('show hbond target a\n')
outcmx.write('cartoon suppress false\nlabel %s  text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\nlabel ontop false\n' %cmxn)
outcmx.write('ui tool show "Side View"\n')

outpml.write("hide labels\n")
if mcount > 2:
	outpml.write('split_states ' + pdbname + '\n')
	for y in range(2,21,1):
		outcmx.write('match #1.%s to #1.1\n' %str(y))
		outpml.write('align %s_%04d, %s_0001\n' %(pdbname,y, pdbname))
outpml.close()
outcmx.close()

