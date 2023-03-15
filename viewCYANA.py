### Mary Clay
import os
import sys
import re
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB',
'CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB',
'CYSSHA':'CA','CYSSHB2':'CB','CYSSHB3':'CB','CYSSQB':'CB',
'ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB',
'GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG',
'PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2',
'GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA',
'HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1',
'HISTHA':'CA','HISTHB2':'CB','HISTHB3':'CB','HISTQB':'CB','HISTHD1':'ND1','HISTHE2':'NE2','HISTHD2':'CD2','HISTHE1':'CE1',
'HIS+HA':'CA','HIS+HB2':'CB','HIS+HB3':'CB','HIS+QB':'CB','HIS+HD1':'ND1','HIS+HE2':'NE2','HIS+HD2':'CD2','HIS+HE1':'CE1',
'ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEQD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1',
'LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ',
'LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUQQD':'CD2,CD1',
'METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHG2':'CG','METHG3':'CG','METQG':'CG','METQE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE',
'ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD2':'ND2',
'PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD',
'CPROHA':'CA','CPROHB2':'CB','CPROHB3':'CB','CPROQB':'CB','CPROHG2':'CG','CPROHG3':'CG','CPROQG':'CG','CPROHD2':'CD','CPROHD3':'CD','CPROQD':'CD',
'GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNQE2':'NE2',
'ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2',
'SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG',
'SEPHA':'CA','SEPHB2':'CB','SEPHB3':'CB','SEPHG':'OG',
'THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2',
'TPOHA':'CA','TPOHB':'CB','TPOHG1':'OG1','TPOHG21':'CG2','TPOHG22':'CG2','TPOHG23':'CG2','TPOQG2':'CG2',
'VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALQQG':'CG1,CG2',
'TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2',
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH',
'PTRHA':'CA','PTRHB2':'CB','PTRHB3':'CB','PTRQB':'CB','PTRQD':'CD1,CD2','PTRQE':'CE1,CE2','PTRHD1':'CD1','PTRHE1':'CE1','PTRHE2':'CE2','PTRHD2':'CD2','PTRHH':'OH',
'ALAH':'N','CYSH':'N','CYSSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','HISTH':'N','HIS+H':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','SEPH':'N','THRH':'N','TPOH':'N','VALH':'N','TRPH':'N','TYRH':'N','PTRH':'N'}

ConTypeDict = {'ALAH':'N', 'ALAHA':'Ali', 'ALAHB':'Methyl', 'ALAHB1':'Methyl', 'ALAHB2':'Methyl', 'ALAHB3':'Methyl', 'ALAQB':'Methyl', 'ALAC':'Other', 'ALACA':'Ali', 'ALACB':'Methyl', 'ALAN':'N', 'CYSSH':'N', 'CYSSHA':'Ali', 'CYSSHB2':'Ali', 'CYSSHB3':'Ali', 'CYSSQB':'Ali', 'CYSSC':'Other', 'CYSSCA':'Ali', 'CYSSCB':'Ali', 'CYSSN':'N', 'CYSH':'N', 'CYSHA':'Ali', 'CYSHB2':'Ali', 'CYSHB3':'Ali', 'CYSQB':'Ali', 'CYSHG':'Ali', 'CYSC':'Other', 'CYSCA':'Ali', 'CYSCB':'Ali', 'CYSN':'N', 'ASPH':'N', 'ASPHA':'Ali', 'ASPHB2':'Ali', 'ASPHB3':'Ali', 'ASPQB':'Ali', 'ASPHD2':'Other', 'ASPC':'Other', 'ASPCA':'Ali', 'ASPCB':'Ali', 'ASPCG':'Ali', 'ASPN':'N', 'GLUH':'N', 'GLUHA':'Ali', 'GLUHB2':'Ali', 'GLUHB3':'Ali', 'GLUQB':'Ali', 'GLUHE2':'Other', 'GLUHG2':'Ali', 'GLUHG3':'Ali', 'GLUQG':'Ali', 'GLUC':'Other', 'GLUCA':'Ali', 'GLUCB':'Ali', 'GLUCG':'Ali', 'GLUCD':'Other', 'GLUN':'N', 'PHEH':'N', 'PHEHA':'Aro', 'PHEHB2':'Aro', 'PHEHB3':'Aro', 'PHEQB':'Aro', 'PHEHD1':'Aro', 'PHEHD2':'Aro', 'PHEQD':'Aro', 'PHEHE1':'Aro', 'PHEHE2':'Aro', 'PHEQE':'Aro', 'PHEHZ':'Aro', 'PHEC':'Other', 'PHECA':'Aro', 'PHECB':'Aro', 'PHECD1':'Aro', 'PHECD2':'Aro', 'PHECE1':'Aro', 'PHECE2':'Aro', 'PHECG':'Aro', 'PHECZ':'Aro', 'PHEN':'N',  'GLYH':'N', 'GLYHA2':'Ali', 'GLYHA3':'Ali', 'GLYC':'Other', 'GLYCA':'Ali', 'GLYN':'N',  'HISH':'N', 'HISHA':'Ali', 'HISHB2':'Ali', 'HISHB3':'Ali', 'HISQB':'Ali', 'HISHD1':'N', 'HISHD2':'Aro', 'HISHE1':'Aro', 'HISC':'Other', 'HISCA':'Ali', 'HISCB':'Ali', 'HISCD2':'Aro', 'HISCE1':'Aro', 'HISCG':'Aro', 'HISN':'N', 'HISND1':'N', 'HISNE2':'N', 'HISTH':'N', 'HISTHA':'Ali', 'HISTHB2':'Ali', 'HISTHB3':'Ali', 'HISTQB':'Ali', 'HISTHD2':'Aro', 'HISTHE1':'Aro', 'HISTHE2':'N', 'HISTC':'Other', 'HISTCA':'Ali', 'HISTCB':'Ali', 'HISTCD2':'Aro', 'HISTCE1':'Aro', 'HISTCG':'Aro', 'HISTN':'N', 'HISTND1':'N', 'HISTNE2':'N', 'HIS+H':'N', 'HIS+HA':'Ali', 'HIS+HB2':'Ali', 'HIS+HB3':'Ali', 'HIS+QB':'Ali', 'HIS+HD1':'N', 'HIS+HD2':'Aro', 'HIS+HE1':'Aro', 'HIS+HE2':'N', 'HIS+C':'Other', 'HIS+CA':'Ali', 'HIS+CB':'Ali', 'HIS+CD2':'Aro', 'HIS+CE1':'Aro', 'HIS+CG':'Aro', 'HIS+N':'N', 'HIS+ND1':'N', 'HIS+NE2':'N', 'ILEH':'N', 'ILEHA':'Ali', 'ILEHB':'Ali', 'ILEHG12':'Ali', 'ILEHG13':'Ali', 'ILEQG1':'Ali', 'ILEHD1':'Methyl', 'ILEHD11':'Methyl', 'ILEHD12':'Methyl', 'ILEHD13':'Methyl', 'ILEQD1':'Methyl', 'ILEHG2':'Methyl', 'ILEHG21':'Methyl', 'ILEHG22':'Methyl', 'ILEHG23':'Methyl', 'ILEQG2':'Methyl', 'ILEC':'Other', 'ILECA':'Ali', 'ILECB':'Ali', 'ILECD1':'Methyl', 'ILECG1':'Ali', 'ILECG2':'Methyl', 'ILEN':'N', 'LYSH':'N', 'LYSHA':'Ali', 'LYSHB2':'Ali', 'LYSHB3':'Ali', 'LYSQB':'Ali', 'LYSHD2':'Ali', 'LYSHD3':'Ali', 'LYSQD':'Ali', 'LYSHE2':'Ali', 'LYSHE3':'Ali', 'LYSQE':'Ali', 'LYSHG2':'Ali', 'LYSHG3':'Ali', 'LYSQG':'Ali', 'LYSC':'Other', 'LYSCA':'Ali', 'LYSCB':'Ali', 'LYSCD':'Ali', 'LYSCE':'Ali', 'LYSCG':'Ali', 'LYSN':'N', 'LYSNZ':'N', 'LYSQZ':'N', 'LYSHZ':'N', 'LEUH':'N', 'LEUHA':'Ali', 'LEUHB2':'Ali', 'LEUHB3':'Ali', 'LEUQB':'Ali', 'LEUHG':'Ali', 'LEUHD1':'Methyl', 'LEUHD11':'Methyl', 'LEUHD12':'Methyl', 'LEUHD13':'Methyl', 'LEUQD1':'Methyl', 'LEUHD2':'Methyl', 'LEUHD21':'Methyl', 'LEUHD22':'Methyl', 'LEUHD23':'Methyl', 'LEUQD2':'Methyl', 'LEUC':'Other', 'LEUCA':'Ali', 'LEUCB':'Ali', 'LEUCG':'Ali', 'LEUCD1':'Methyl', 'LEUCD2':'Methyl', 'LEUN':'N', 'METH':'N', 'METHA':'Ali', 'METHB2':'Ali', 'METHB3':'Ali', 'METQB':'Ali', 'METHG2':'Ali', 'METHG3':'Ali', 'METQG':'Ali', 'METHE':'Methyl', 'METHE1':'Methyl', 'METHE2':'Methyl', 'METHE3':'Methyl', 'METQE':'Methyl', 'METC':'Other', 'METCA':'Ali', 'METCB':'Ali', 'METCE':'Methyl', 'METCG':'Ali', 'METN':'N', 'ASNH':'N', 'ASNHA':'Ali', 'ASNHB2':'Ali', 'ASNHB3':'Ali', 'ASNQB':'Ali', 'ASNHD21':'N', 'ASNHD22':'N', 'ASNQD':'N', 'ASNC':'Other', 'ASNCA':'Ali', 'ASNCB':'Ali', 'ASNCG':'Other', 'ASNN':'N', 'ASNND2':'N', 'PROHA':'Ali', 'PROHB2':'Ali', 'PROHB3':'Ali', 'PROQB':'Ali', 'PROHD2':'Ali', 'PROHD3':'Ali', 'PROQD':'Ali', 'PROHG2':'Ali', 'PROHG3':'Ali', 'PROQG':'Ali', 'PROC':'Other', 'PROCA':'Ali', 'PROCB':'Ali', 'PROCD':'Ali', 'PROCG':'Ali', 'PRON':'N', 'CPROHA':'Ali', 'CPROHB2':'Ali', 'CPROHB3':'Ali', 'CPROQB':'Ali', 'CPROHD2':'Ali', 'CPROHD3':'Ali', 'CPROQD':'Ali', 'CPROHG2':'Ali', 'CPROHG3':'Ali', 'CPROQG':'Ali', 'CPROC':'Other', 'CPROCA':'Ali', 'CPROCB':'Ali', 'CPROCD':'Ali', 'CPROCG':'Ali', 'CPRON':'N', 'GLNH':'N', 'GLNHA':'Ali', 'GLNHB2':'Ali', 'GLNHB3':'Ali', 'GLNQB':'Ali', 'GLNHE21':'N', 'GLNHE22':'N', 'GLNQE2':'N', 'GLNHG2':'Ali', 'GLNHG3':'Ali', 'GLNQG':'Ali', 'GLNC':'Other', 'GLNCA':'Ali', 'GLNCB':'Ali', 'GLNCD':'Other', 'GLNCG':'Ali', 'GLNN':'N', 'GLNNE2':'N', 'ARGH':'N', 'ARGHA':'Ali', 'ARGHB2':'Ali', 'ARGHB3':'Ali', 'ARGQB':'Ali', 'ARGHD2':'Ali', 'ARGHD3':'Ali', 'ARGQD':'Ali', 'ARGHG2':'Ali', 'ARGHG3':'Ali', 'ARGQG':'Ali', 'ARGHH11':'N', 'ARGHH12':'N', 'ARGQH1':'N', 'ARGHH21':'N', 'ARGHH22':'N', 'ARGQH2':'N', 'ARGC':'Other', 'ARGCA':'Ali', 'ARGCB':'Ali', 'ARGCD':'Ali', 'ARGCG':'Ali', 'ARGCZ':'Ali', 'ARGN':'N', 'ARGNE':'N', 'ARGNH1':'N', 'ARGNH2':'N', 'ARGHE':'N', 'SERH':'N', 'SERHA':'Ali', 'SERHB2':'Ali', 'SERHB3':'Ali', 'SERQB':'Ali', 'SERHG':'Other', 'SERC':'Other', 'SERCA':'Ali', 'SERCB':'Ali', 'SERN':'N', 'SEPH':'N', 'SEPHA':'Ali', 'SEPHB2':'Ali', 'SEPHB3':'Ali', 'SEPQB':'Ali', 'SEPHG':'Other', 'SEPC':'Other', 'SEPCA':'Ali', 'SEPCB':'Ali', 'SEPN':'N', 'THRH':'N', 'THRHA':'Ali', 'THRHB':'Ali', 'THRHG1':'Other', 'THRHG2':'Methyl', 'THRQG2':'Methyl', 'THRC':'Other', 'THRCA':'Ali', 'THRCB':'Ali', 'THRCG2':'Methyl', 'THRN':'N', 'TPOH':'N', 'TPOHA':'Ali', 'TPOHB':'Ali', 'TPOHG1':'Other', 'TPOHG2':'Methyl', 'TPOQG2':'Methyl', 'TPOC':'Other', 'TPOCA':'Ali', 'TPOCB':'Ali', 'TPOCG2':'Methyl', 'TPON':'N', 'VALH':'N', 'VALHA':'Ali', 'VALHB':'Ali', 'VALHG1':'Methyl', 'VALQG1':'Methyl', 'VALHG2':'Methyl', 'VALQG2':'Methyl', 'VALC':'Other', 'VALCA':'Ali', 'VALCB':'Ali', 'VALCG1':'Methyl', 'VALCG2':'Methyl', 'VALN':'N', 'TRPH':'N', 'TRPHA':'Ali', 'TRPHB2':'Ali', 'TRPHB3':'Ali', 'TRPQB':'Ali', 'TRPHD1':'Aro', 'TRPHE1':'N', 'TRPHE3':'Aro', 'TRPHH2':'Aro', 'TRPHZ2':'Aro', 'TRPHZ3':'Aro', 'TRPC':'Other', 'TRPCA':'Ali', 'TRPCB':'Ali', 'TRPCD1':'Aro', 'TRPCD2':'Aro', 'TRPCE2':'Aro', 'TRPCE3':'Aro', 'TRPCG':'Aro', 'TRPCH2':'Aro', 'TRPCZ2':'Aro', 'TRPCZ3':'Aro', 'TRPN':'N', 'TRPNE1':'N', 'TYRH':'N', 'TYRHA':'Aro', 'TYRHB2':'Aro', 'TYRHB3':'Aro', 'TYRQB':'Aro', 'TYRHD1':'Aro', 'TYRHD2':'Aro', 'TYRQD':'Aro', 'TYRHE1':'Aro', 'TYRHE2':'Aro', 'TYRQE':'Aro', 'TYRHH':'Other', 'TYRC':'Other', 'TYRCA':'Aro', 'TYRCB':'Aro', 'TYRCD1':'Aro', 'TYRCD2':'Aro', 'TYRCE1':'Aro', 'TYRCE2':'Aro', 'TYRCG':'Aro', 'TYRCZ':'Aro', 'TYRN':'N', 'PTRH':'N', 'PTRHA':'Aro', 'PTRHB2':'Aro', 'PTRHB3':'Aro', 'PTRQB':'Aro', 'PTRHD1':'Aro', 'PTRHD2':'Aro', 'PTRQD':'Aro', 'PTRHE1':'Aro', 'PTRHE2':'Aro', 'PTRQE':'Aro', 'PTRHH':'Other', 'PTRC':'Other', 'PTRCA':'Aro', 'PTRCB':'Aro', 'PTRCD1':'Aro', 'PTRCD2':'Aro', 'PTRCE1':'Aro', 'PTRCE2':'Aro', 'PTRCG':'Aro', 'PTRCZ':'Aro', 'PTRN':'N', 'N-N':'N_N', 'Aro-Aro':'Aro_Aro', 'Methyl-Methyl':'Methyl_Methyl', 'Ali-Ali':'Ali_Ali', 'Aro-Methyl':'Methyl_Aro', 'Methyl-Aro':'Methyl_Aro', 'Methyl-N':'N_Methyl', 'N-Methyl':'N_Methyl', 'Aro-N':'N_Aro', 'N-Aro':'N_Aro', 'Methyl-Ali':'Ali_Methyl', 'Ali-Methyl':'Ali_Methyl', 'Aro-Ali':'Ali_Aro', 'Ali-Aro':'Ali_Aro', 'Ali-N':'N_Ali', 'N-Ali':'N_Ali', 'Other-Other':'Other','Other-N':'Other','N-Other':'Other','Other-Ali':'Other','Ali-Other':'Other','Other-Methyl':'Other','Methyl-Other':'Other', 'Other-Aro':'Other', 'Aro-Other':'Other'}

# ['ALAH','CYSH','ASPH','GLUH','PHEH','GLYH','HISH','ILEH','LYSH','LEUH','METH','ASNH','GLNH','ARGH','SERH','THRH','VALH','TRPH','TYRH']
if len(sys.argv)==1:
	print('''

Usage: 
	viewcya [pdb] [TALOS]

	viewcya AF-FGFR3_KD.pdb ../TALOS2

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
colors = ['royalblue','forest','yellowgreen', 'darkorange','purple','lightseagreen ','darkkhaki','peru','saddlebrown','mediumpurple','blue']
ConectionTypes = ['N_N','N_Methyl','N_Aro','Methyl_Methyl','Methyl_Aro','Aro_Aro','N_Ali','Ali_Ali','Ali_Aro','Ali_Methyl','Other']

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
if not os.path.exists(outdir +'pseudobonds/'):
	os.makedirs(outdir +'pseudobonds/')

## open the CALC.cya file to get the peaks list and additonal constraint files used in the calculation. 
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
lols = [con for con in manualongcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualongcons if 'aco' in con]
mn = 1
outpml = open(outdir + 'CYANA_input.pml','w')
outpml.write('load '+pdb_path+'\n')
outpml.write('load '+pymol_pdb_path+'\n')
outpml.write('set_color royalblue = [65,105,225]\nset_color forest = [34,139,34]\nset_color yellowgreen = [154,205,50]\nset_color darkorange = [255,140,0]\nset_color purple = [128,0,128]\nset_color lightseagreen = [32,178,170]\nset_color darkkhaki = [189,183,107]\nset_color peru = [205,133,63]\nset_color saddlebrown = [139,69,19]\nset_color gold = [255,215,0]\nset_color navy = [0,0,128]\nset_color darkturquoise = [0,206,209]\nset_color pink = [255,192,203]\nset_color cyan = [0,255,255]\nset_color paleturquoise = [175,238,238]\nset_color lightsalmon = [255,160,122]\nset_color khaki = [240,230,140]\nset_color yellowgreen = [154,205,50]\nset_color thistle = [216,191,216]\nset_color aquamarine = [127,255,212]\nset_color plum = [221,160,221]\nset_color lightpink = [255,182,193]\nset_color mediumvioletred = [199,21,133]\nset_color firebrick = [178,34,34]\nset_color lightcoral = [240,128,128]\nset_color deeppink = [255,20,147]\nset_color hotpink = [255,105,180]\nset_color mediumpurple = [147,112,219]\nset_color navy = [0,0,128]\nset_color cornflowerblue = [100,149,237]\n')
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
hbond = open(outdir + 'pseudobonds/hbond.pb','w')
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
mn+=2
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
mn+=1
for uplfile in upls:
	for conect in ConectionTypes:
		exec("{:}_{:}_pb = []".format(uplfile.replace('.upl',''),conect))
		exec("group{:}{:} = 'group {:}_{:}, '".format(uplfile.replace('.upl',''),conect, uplfile.replace('.upl',''), conect))
for uplfile in upls:
	fin = open(uplfile,'r')
	# outpb.write("; halfbond = false\n; color = {:}\n; radius = 0.15\n; dashes = 10\n".format(colors[x]))
	for line in fin.readlines():
		if line.split():
			cns = line.split()
			if "#" not in cns[0]:
				if cns[1]+cns[2] in ConTypeDict.keys():ct1 = ConTypeDict[cns[1]+cns[2]]
				if cns[1]+cns[2] not in ConTypeDict.keys(): ct1 = 'Other'
				if cns[4]+cns[5] in ConTypeDict.keys():ct2 = ConTypeDict[cns[4]+cns[5]]
				if cns[4]+cns[5] not in ConTypeDict.keys(): ct2 = 'Other'
				ctype = ConTypeDict["{:}-{:}".format(ct1,ct2)]
				outpb = eval('{:}_{:}_pb'.format(uplfile.replace('.upl',''), ctype))
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
						outpb.append('{:}:{:}@{:} {:}:{:}@{:}\n'.format(cmxn, cns[0], atom1, cmxn, cns[3],atom2))
						outpml.write('distance UPL{:}, {:} and resi {:} and name {:}, {:} and resi {:} and name {:}\n'.format(str(u), pmln, cns[0], atom1, pmln, cns[3], atom2))
						exec('group{:}{:} = group{:}{:} + "UPL{:}"'.format(uplfile.replace('.upl',''), ctype,uplfile.replace('.upl',''), ctype,u))
				if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
					sidelist.append(cns[0])
				if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
					sidelist.append(cns[3])

for uplfile in upls:
	for x in range(len(ConectionTypes)):
		pbs = eval('{:}_{:}_pb'.format(uplfile.replace('.upl',''),ConectionTypes[x]))
		if len(pbs) > 1:
			mn+=1
			pbout = open('{:}pseudobonds/{:}_{:}.pb'.format(outdir, uplfile.replace('.upl',''), ConectionTypes[x]),'w')
			pbout.write("; halfbond = false\n; color = " + colors[x] + "\n; radius = 0.1\n; dashes = 0\n")
			pbout.writelines(pbs)
			outcmx.write('open pseudobonds/{:}_{:}.pb\n'.format(uplfile.replace('.upl',''), ConectionTypes[x]))
			outcmx.write('color #{:} {:}\n'.format(str(mn),colors[x]))
			groupstr = eval('group{:}{:}'.format(uplfile.replace('.upl',''),ConectionTypes[x]))
			outpml.write(groupstr + '\n')
			outpml.write('color {:}, {:}_{:}\n'.format(colors[x],uplfile.replace('.upl',''),ConectionTypes[x]))


sidechains = 'show #2:'
for res in sidelist:
	sidechains = sidechains + res + ','
outcmx.write(sidechains[:-1] + '\n')
outpml.write('show sticks, {:} and resi {:}\n'.format(pdbname, sidechains[sidechains.index(':')+1:-1].replace(',','+')))

outpml.write('color pink, hbond\n')
selhbond = selhbond[:-1] + '@O,N\n'
sidehbond = sidehbond[:-1] + '\n'
outcmx.write('open pseudobonds/hbond.pb\n')
mn+=1
outcmx.write('color #{:} {:}\n'.format(str(mn),'pink'))
outcmx.write(selhbond)
outcmx.write('show hbond\n')
if ',' in sidehbond:
	outcmx.write(sidehbond)
	outcmx.write('show shbond  target a\nhide H\nshow hbond target a\n')


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
outcmx.write('hide H\nshow {:}@N,H target a\n'.format(cmxn))

outcmx.write('hide #{:}@C,O,N,H\n'.format(3))

outcmx.write('ui tool show "Side View"\n')
outpml.write("hide labels\n")
for y in range(2,21,1):
	outpml.write('delete {:}_{:04d}\n'.format(pdbname,y))
outpml.close()
outcmx.close()


