### Mary Clay
import os
import sys
import re
import pandas as pd
import numpy as np
import GetDihe as Dihed
import noa_analysis as noaa
from datetime import datetime
import glob

# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y-%m-%d %H:%M")

Vnum = '2.4'
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB','CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB','CYSSHA':'CA','CYSSHB2':'CB','CYSSHB3':'CB','CYSSQB':'CB','ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB','GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG','PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2','GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA','HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1','HISTHA':'CA','HISTHB2':'CB','HISTHB3':'CB','HISTQB':'CB','HISTHD1':'ND1','HISTHE2':'NE2','HISTHD2':'CD2','HISTHE1':'CE1','HIS+HA':'CA','HIS+HB2':'CB','HIS+HB3':'CB','HIS+QB':'CB','HIS+HD1':'ND1','HIS+HE2':'NE2','HIS+HD2':'CD2','HIS+HE1':'CE1','ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEQD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1','LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ','LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUQQD':'CD2,CD1','METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHG2':'CG','METHG3':'CG','METQG':'CG','METQE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE','ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD21':'ND2','ASNQD22':'ND2','ASNQD2':'ND2','PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD','CPROHA':'CA','CPROHB2':'CB','CPROHB3':'CB','CPROQB':'CB','CPROHG2':'CG','CPROHG3':'CG','CPROQG':'CG','CPROHD2':'CD','CPROHD3':'CD','CPROQD':'CD','GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNQE21':'NE2','GLNQE22':'NE2','GLNQE2':'NE2','ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2','SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG','SEPHA':'CA','SEPHB2':'CB','SEPHB3':'CB','SEPHG':'OG','THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2','TPOHA':'CA','TPOHB':'CB','TPOHG1':'OG1','TPOHG21':'CG2','TPOHG22':'CG2','TPOHG23':'CG2','TPOQG2':'CG2','VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALQQG':'CG1,CG2','TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2','TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH','PTRHA':'CA','PTRHB2':'CB','PTRHB3':'CB','PTRQB':'CB','PTRQD':'CD1,CD2','PTRQE':'CE1,CE2','PTRHD1':'CD1','PTRHE1':'CE1','PTRHE2':'CE2','PTRHD2':'CD2','PTRHH':'OH','ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N','PTRH':'N','HISTH':'N','HISEH':'N','TPOH':'N','SEPH':'N'}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N','PTRH':'N','HISTH':'N','HISEH':'N','TPOH':'N','SEPH':'N'}
ConTypeDict = {'ALAH':'N', 'ALAHA':'Ali', 'ALAHB':'Methyl', 'ALAHB1':'Methyl', 'ALAHB2':'Methyl', 'ALAHB3':'Methyl', 'ALAQB':'Methyl', 'ALAC':'Other', 'ALACA':'Ali', 'ALACB':'Methyl', 'ALAN':'N', 'CYSSH':'N', 'CYSSHA':'Ali', 'CYSSHB2':'Ali', 'CYSSHB3':'Ali', 'CYSSQB':'Ali', 'CYSSC':'Other', 'CYSSCA':'Ali', 'CYSSCB':'Ali', 'CYSSN':'N', 'CYSH':'N', 'CYSHA':'Ali', 'CYSHB2':'Ali', 'CYSHB3':'Ali', 'CYSQB':'Ali', 'CYSHG':'Ali', 'CYSC':'Other', 'CYSCA':'Ali', 'CYSCB':'Ali', 'CYSN':'N', 'ASPH':'N', 'ASPHA':'Ali', 'ASPHB2':'Ali', 'ASPHB3':'Ali', 'ASPQB':'Ali', 'ASPHD2':'Other', 'ASPC':'Other', 'ASPCA':'Ali', 'ASPCB':'Ali', 'ASPCG':'Ali', 'ASPN':'N', 'GLUH':'N', 'GLUHA':'Ali', 'GLUHB2':'Ali', 'GLUHB3':'Ali', 'GLUQB':'Ali', 'GLUHE2':'Other', 'GLUHG2':'Ali', 'GLUHG3':'Ali', 'GLUQG':'Ali', 'GLUC':'Other', 'GLUCA':'Ali', 'GLUCB':'Ali', 'GLUCG':'Ali', 'GLUCD':'Other', 'GLUN':'N', 'PHEH':'N', 'PHEHA':'Aro', 'PHEHB2':'Aro', 'PHEHB3':'Aro', 'PHEQB':'Aro', 'PHEHD1':'Aro', 'PHEHD2':'Aro', 'PHEQD':'Aro', 'PHEHE1':'Aro', 'PHEHE2':'Aro', 'PHEQE':'Aro', 'PHEHZ':'Aro', 'PHEC':'Other', 'PHECA':'Aro', 'PHECB':'Aro', 'PHECD1':'Aro', 'PHECD2':'Aro', 'PHECE1':'Aro', 'PHECE2':'Aro', 'PHECG':'Aro', 'PHECZ':'Aro', 'PHEN':'N',  'GLYH':'N', 'GLYHA2':'Ali', 'GLYHA3':'Ali', 'GLYC':'Other', 'GLYCA':'Ali', 'GLYN':'N',  'HISH':'N', 'HISHA':'Ali', 'HISHB2':'Ali', 'HISHB3':'Ali', 'HISQB':'Ali', 'HISHD1':'N', 'HISHD2':'Aro', 'HISHE1':'Aro', 'HISC':'Other', 'HISCA':'Ali', 'HISCB':'Ali', 'HISCD2':'Aro', 'HISCE1':'Aro', 'HISCG':'Aro', 'HISN':'N', 'HISND1':'N', 'HISNE2':'N', 'HISTH':'N', 'HISTHA':'Ali', 'HISTHB2':'Ali', 'HISTHB3':'Ali', 'HISTQB':'Ali', 'HISTHD2':'Aro', 'HISTHE1':'Aro', 'HISTHE2':'N', 'HISTC':'Other', 'HISTCA':'Ali', 'HISTCB':'Ali', 'HISTCD2':'Aro', 'HISTCE1':'Aro', 'HISTCG':'Aro', 'HISTN':'N', 'HISTND1':'N', 'HISTNE2':'N', 'HIS+H':'N', 'HIS+HA':'Ali', 'HIS+HB2':'Ali', 'HIS+HB3':'Ali', 'HIS+QB':'Ali', 'HIS+HD1':'N', 'HIS+HD2':'Aro', 'HIS+HE1':'Aro', 'HIS+HE2':'N', 'HIS+C':'Other', 'HIS+CA':'Ali', 'HIS+CB':'Ali', 'HIS+CD2':'Aro', 'HIS+CE1':'Aro', 'HIS+CG':'Aro', 'HIS+N':'N', 'HIS+ND1':'N', 'HIS+NE2':'N', 'ILEH':'N', 'ILEHA':'Ali', 'ILEHB':'Ali', 'ILEHG12':'Ali', 'ILEHG13':'Ali', 'ILEQG1':'Ali', 'ILEHD1':'Methyl', 'ILEHD11':'Methyl', 'ILEHD12':'Methyl', 'ILEHD13':'Methyl', 'ILEQD1':'Methyl', 'ILEHG2':'Methyl', 'ILEHG21':'Methyl', 'ILEHG22':'Methyl', 'ILEHG23':'Methyl', 'ILEQG2':'Methyl', 'ILEC':'Other', 'ILECA':'Ali', 'ILECB':'Ali', 'ILECD1':'Methyl', 'ILECG1':'Ali', 'ILECG2':'Methyl', 'ILEN':'N', 'LYSH':'N', 'LYSHA':'Ali', 'LYSHB2':'Ali', 'LYSHB3':'Ali', 'LYSQB':'Ali', 'LYSHD2':'Ali', 'LYSHD3':'Ali', 'LYSQD':'Ali', 'LYSHE2':'Ali', 'LYSHE3':'Ali', 'LYSQE':'Ali', 'LYSHG2':'Ali', 'LYSHG3':'Ali', 'LYSQG':'Ali', 'LYSC':'Other', 'LYSCA':'Ali', 'LYSCB':'Ali', 'LYSCD':'Ali', 'LYSCE':'Ali', 'LYSCG':'Ali', 'LYSN':'N', 'LYSNZ':'N', 'LYSQZ':'N', 'LYSHZ':'N', 'LEUH':'N', 'LEUHA':'Ali', 'LEUHB2':'Ali', 'LEUHB3':'Ali', 'LEUQB':'Ali', 'LEUHG':'Ali', 'LEUHD1':'Methyl', 'LEUHD11':'Methyl', 'LEUHD12':'Methyl', 'LEUHD13':'Methyl', 'LEUQD1':'Methyl', 'LEUHD2':'Methyl', 'LEUHD21':'Methyl', 'LEUHD22':'Methyl', 'LEUHD23':'Methyl', 'LEUQD2':'Methyl', 'LEUC':'Other', 'LEUCA':'Ali', 'LEUCB':'Ali', 'LEUCG':'Ali', 'LEUCD1':'Methyl', 'LEUCD2':'Methyl', 'LEUN':'N', 'METH':'N', 'METHA':'Ali', 'METHB2':'Ali', 'METHB3':'Ali', 'METQB':'Ali', 'METHG2':'Ali', 'METHG3':'Ali', 'METQG':'Ali', 'METHE':'Methyl', 'METHE1':'Methyl', 'METHE2':'Methyl', 'METHE3':'Methyl', 'METQE':'Methyl', 'METC':'Other', 'METCA':'Ali', 'METCB':'Ali', 'METCE':'Methyl', 'METCG':'Ali', 'METN':'N', 'ASNH':'N', 'ASNHA':'Ali', 'ASNHB2':'Ali', 'ASNHB3':'Ali', 'ASNQB':'Ali', 'ASNHD21':'N', 'ASNHD22':'N', 'ASNQD':'N', 'ASNC':'Other', 'ASNCA':'Ali', 'ASNCB':'Ali', 'ASNCG':'Other', 'ASNN':'N', 'ASNND2':'N', 'PROHA':'Ali', 'PROHB2':'Ali', 'PROHB3':'Ali', 'PROQB':'Ali', 'PROHD2':'Ali', 'PROHD3':'Ali', 'PROQD':'Ali', 'PROHG2':'Ali', 'PROHG3':'Ali', 'PROQG':'Ali', 'PROC':'Other', 'PROCA':'Ali', 'PROCB':'Ali', 'PROCD':'Ali', 'PROCG':'Ali', 'PRON':'N', 'CPROHA':'Ali', 'CPROHB2':'Ali', 'CPROHB3':'Ali', 'CPROQB':'Ali', 'CPROHD2':'Ali', 'CPROHD3':'Ali', 'CPROQD':'Ali', 'CPROHG2':'Ali', 'CPROHG3':'Ali', 'CPROQG':'Ali', 'CPROC':'Other', 'CPROCA':'Ali', 'CPROCB':'Ali', 'CPROCD':'Ali', 'CPROCG':'Ali', 'CPRON':'N', 'GLNH':'N', 'GLNHA':'Ali', 'GLNHB2':'Ali', 'GLNHB3':'Ali', 'GLNQB':'Ali', 'GLNHE21':'N', 'GLNHE22':'N', 'GLNQE2':'N', 'GLNHG2':'Ali', 'GLNHG3':'Ali', 'GLNQG':'Ali', 'GLNC':'Other', 'GLNCA':'Ali', 'GLNCB':'Ali', 'GLNCD':'Other', 'GLNCG':'Ali', 'GLNN':'N', 'GLNNE2':'N', 'ARGH':'N', 'ARGHA':'Ali', 'ARGHB2':'Ali', 'ARGHB3':'Ali', 'ARGQB':'Ali', 'ARGHD2':'Ali', 'ARGHD3':'Ali', 'ARGQD':'Ali', 'ARGHG2':'Ali', 'ARGHG3':'Ali', 'ARGQG':'Ali', 'ARGHH11':'N', 'ARGHH12':'N', 'ARGQH1':'N', 'ARGHH21':'N', 'ARGHH22':'N', 'ARGQH2':'N', 'ARGC':'Other', 'ARGCA':'Ali', 'ARGCB':'Ali', 'ARGCD':'Ali', 'ARGCG':'Ali', 'ARGCZ':'Ali', 'ARGN':'N', 'ARGNE':'N', 'ARGNH1':'N', 'ARGNH2':'N', 'ARGHE':'N', 'SERH':'N', 'SERHA':'Ali', 'SERHB2':'Ali', 'SERHB3':'Ali', 'SERQB':'Ali', 'SERHG':'Other', 'SERC':'Other', 'SERCA':'Ali', 'SERCB':'Ali', 'SERN':'N', 'SEPH':'N', 'SEPHA':'Ali', 'SEPHB2':'Ali', 'SEPHB3':'Ali', 'SEPQB':'Ali', 'SEPHG':'Other', 'SEPC':'Other', 'SEPCA':'Ali', 'SEPCB':'Ali', 'SEPN':'N', 'THRH':'N', 'THRHA':'Ali', 'THRHB':'Ali', 'THRHG1':'Other', 'THRHG2':'Methyl', 'THRQG2':'Methyl', 'THRC':'Other', 'THRCA':'Ali', 'THRCB':'Ali', 'THRCG2':'Methyl', 'THRN':'N', 'TPOH':'N', 'TPOHA':'Ali', 'TPOHB':'Ali', 'TPOHG1':'Other', 'TPOHG2':'Methyl', 'TPOQG2':'Methyl', 'TPOC':'Other', 'TPOCA':'Ali', 'TPOCB':'Ali', 'TPOCG2':'Methyl', 'TPON':'N', 'VALH':'N', 'VALHA':'Ali', 'VALHB':'Ali', 'VALHG1':'Methyl', 'VALQG1':'Methyl', 'VALHG2':'Methyl', 'VALQG2':'Methyl', 'VALC':'Other', 'VALCA':'Ali', 'VALCB':'Ali', 'VALCG1':'Methyl', 'VALCG2':'Methyl', 'VALN':'N', 'TRPH':'N', 'TRPHA':'Ali', 'TRPHB2':'Ali', 'TRPHB3':'Ali', 'TRPQB':'Ali', 'TRPHD1':'Aro', 'TRPHE1':'N', 'TRPHE3':'Aro', 'TRPHH2':'Aro', 'TRPHZ2':'Aro', 'TRPHZ3':'Aro', 'TRPC':'Other', 'TRPCA':'Ali', 'TRPCB':'Ali', 'TRPCD1':'Aro', 'TRPCD2':'Aro', 'TRPCE2':'Aro', 'TRPCE3':'Aro', 'TRPCG':'Aro', 'TRPCH2':'Aro', 'TRPCZ2':'Aro', 'TRPCZ3':'Aro', 'TRPN':'N', 'TRPNE1':'N', 'TYRH':'N', 'TYRHA':'Aro', 'TYRHB2':'Aro', 'TYRHB3':'Aro', 'TYRQB':'Aro', 'TYRHD1':'Aro', 'TYRHD2':'Aro', 'TYRQD':'Aro', 'TYRHE1':'Aro', 'TYRHE2':'Aro', 'TYRQE':'Aro', 'TYRHH':'Other', 'TYRC':'Other', 'TYRCA':'Aro', 'TYRCB':'Aro', 'TYRCD1':'Aro', 'TYRCD2':'Aro', 'TYRCE1':'Aro', 'TYRCE2':'Aro', 'TYRCG':'Aro', 'TYRCZ':'Aro', 'TYRN':'N', 'PTRH':'N', 'PTRHA':'Aro', 'PTRHB2':'Aro', 'PTRHB3':'Aro', 'PTRQB':'Aro', 'PTRHD1':'Aro', 'PTRHD2':'Aro', 'PTRQD':'Aro', 'PTRHE1':'Aro', 'PTRHE2':'Aro', 'PTRQE':'Aro', 'PTRHH':'Other', 'PTRC':'Other', 'PTRCA':'Aro', 'PTRCB':'Aro', 'PTRCD1':'Aro', 'PTRCD2':'Aro', 'PTRCE1':'Aro', 'PTRCE2':'Aro', 'PTRCG':'Aro', 'PTRCZ':'Aro', 'PTRN':'N', 'N-N':'N_N', 'Aro-Aro':'Aro_Aro', 'Methyl-Methyl':'Methyl_Methyl', 'Ali-Ali':'Ali_Ali', 'Aro-Methyl':'Methyl_Aro', 'Methyl-Aro':'Methyl_Aro', 'Methyl-N':'N_Methyl', 'N-Methyl':'N_Methyl', 'Aro-N':'N_Aro', 'N-Aro':'N_Aro', 'Methyl-Ali':'Ali_Methyl', 'Ali-Methyl':'Ali_Methyl', 'Aro-Ali':'Ali_Aro', 'Ali-Aro':'Ali_Aro', 'Ali-N':'N_Ali', 'N-Ali':'N_Ali', 'Other-Other':'Other','Other-N':'Other','N-Other':'Other','Other-Ali':'Other','Ali-Other':'Other','Other-Methyl':'Other','Methyl-Other':'Other', 'Other-Aro':'Other', 'Aro-Other':'Other'}
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}
Ambiguous = {'ARGQB':'HB2,HB3','ARGQG':'HG2,HG3','ARGQG':'HG2,HG3','ARGQH1':'HH11,HH12','ARGQH2':'HH21,HH22','ASNQB':'HB2,HB3','ASNQD2':'HD21,HD22','ASPQB':'HB2,HB3','CYSQB':'HB2,HB3','CYSSQB':'HB2,HB3','GLNQB':'HB2,HB3','GLNQG':'HG2,HG3','GLNQE2':'HE21,HE22','GLUQB':'HB2,HB3','GLUQG':'HG2,HG3','GLYQA':'HA2,HA3','HISQB':'HB2,HB3','HISTQB':'HB2,HB3','ILEQG1':'HG12,HG13','LEUQB':'HB2,HB3','LEUQQD':'QD1,QD2','LYSQB':'HB2,HB3','LYSQG':'HG2,HG3','LYSQD':'HD2,HD3','LYSQE':'HE2,HE3','METQB':'HB2,HB3','METQG':'HG2,HG3','PHEQB':'HB2,HB3','PHEQG':'HG2,HG3','PHEQD':'HD2,HD3','PHEQE':'HE1,HE2','PROQB':'HB2,HB3','PROQG':'HG2,HG3','PROQD':'HD2,HD3','SERQB':'HB2,HB3','TRPQB':'HB2,HB3','TYRQB':'HB2,HB3','TYRQG':'HG2,HG3','TYRQD':'HD2,HD3','TYRQE':'HE1,HE2','PTRQB':'HB2,HB3','PTRQG':'HG2,HG3','PTRQD':'HD2,HD3','PTRQE':'HE1,HE2','VALQQG':'QG1,QG2'}

if len(sys.argv)==1:
  print('''

Usage: 
  cyanapra [pdb] [upl]

Required Input:

  PDB         PDB to be used typically the final.pdb or pdb after CNS refinement
        If this is not located in current directory provide path
          CNS/refinePDB/r12_cya.pdb

  upl             What upl file would you like to use? final.upl cycle?.upl
          Which ever upl you specify determines the overview file used. 

OutPut:
  post_cyana_analysis/
    pdb.cxc
    pdb.pml
    upl_summary.txt
    upl_overview.pdf
    pdb_Phi.csv, pdb_Psi.csv, pdb_Chi1.csv, pdb_Chi2.csv 
    Annotated upl files 
    pseudobonds/ (Groups of Distance restraints)
      upl_constraint_type.pb
      upl_poor.pb: restraints with SUP < 0.5
      upl_short.pb:  non intramolecular distance restraints d < 3.0
      upl_long.pb: Long distance restraints d > 6.0
      upl_viol_upls.pb: UPLs violated in 10 or more structures 
      input.pb(s): input upl file(s) specified in CALC.cya
      hbond.pb
  noa_analysis/ 
    Assignment_Summary.txt 
    peak-list.list
''')
  exit()
colors = ['royalblue','forest','yellowgreen', 'darkorange','purple','lightseagreen ','darkkhaki','peru','saddlebrown','mediumpurple','blue','purple']
ConectionTypes = ['N_N','N_Methyl','N_Aro','Methyl_Methyl','Methyl_Aro','Aro_Aro','N_Ali','Ali_Ali','Ali_Aro','Ali_Methyl','Other','intermolecular']
## ---------------------------------------------------------------------------
## Read in termina input and setup output directories 
## ---------------------------------------------------------------------------
MasterDict = {}
cwd = os.getcwd() + '/'
outdir = cwd + 'post_cyana_analysis/'
noadir = cwd  +'noa_analysis/'
in_pdb = sys.argv[1]
fupl = sys.argv[2]
pdbname = in_pdb.split('/')[-1].split('.')[0]
fovw = fupl.replace('.upl','.ovw')
calc = cwd + 'CALC.cya'
outname = in_pdb.split('.')[0]
init = cwd + 'init.cya'
seq = [line.strip().split() for line in open(cwd + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0] and line.split()[0] in AAA_dict.keys()]
Seqdict = {}
Sequence, ASequence = [],[]
for resn,resi in seq:
  if resn in AAA_dict.keys():
    Seqdict[resi] = AAA_dict[resn] + resi
    ASequence.append(AAA_dict[resn] + resi)
    Sequence.append(resi)
upldf = pd.DataFrame(index = ASequence, columns=['cya','long','viol','input','found input','viol input'])
upldf['cya'] = np.zeros(len(Sequence))
upldf['long'] = np.zeros(len(Sequence))
upldf['viol'] = np.zeros(len(Sequence))
upldf['input'] = np.zeros(len(Sequence))
upldf['found input'] = np.zeros(len(Sequence))
upldf['viol input'] = np.zeros(len(Sequence))
## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
  os.makedirs(outdir)
if not os.path.exists(outdir +'pseudobonds/'):
  os.makedirs(outdir +'pseudobonds/')
if not os.path.exists(noadir):
  os.makedirs(noadir)
checkcons = open(outdir + outname + '_summary.txt','w')
## open the CALC.cya file to get the peaks list and additional constraint files used in the calculation. 
cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line and not re.match(r'^\s*#', line)][0].split()[2].split(',')
lengths,upllengths = [], []
for plist in cya_plists:
  lengths.append(len(plist))
pad = ''
for x in range(max(lengths)):
  pad = pad + ' '

manualcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualcons if 'upl' in con and 'hbond' not in con]
hbonds = [con for con in manualcons if 'upl' in con and 'hbond' in con]
lols = [con for con in manualcons if 'lol' in con and 'hbond' not in con]
dihed = [con for con in manualcons if 'aco' in con if con != 'inital.aco']
noa = cwd + 'cycle7.noa'
noalines = open(noa).readlines()

shortsum = open(outdir + 'Short_stats.txt','w')
checkcons.write('## Generated using CYANApra_{:} on {:} \n'.format(Vnum,dt_string))
shortsum.write('## Generated using CYANApra_{:} on {:} \n'.format(Vnum,dt_string))
peaklisttxt = '{:}                                         Assignments \n{:}#peaks   upl  Viol Unique  Multiple  Unused  None  Diagonal Increased upl    %\n'.format(pad,pad)
#checkcons.write('{:}                                         Assignments \n{:}#peaks   upl  Viol Unique  Multiple  Unused  None  Diagonal Increased upl    %\n'.format(pad,pad))
#print('{:}                                         Assignments \n{:}#peaks   upl  Viol Unique  Multiple  Unused  None  Diagonal Increased upl    %'.format(pad,pad))
shortsumtext = '{:} #peaks  % Assignment\n'.format(pad)

## ---------------------------------------------------------------------------
## Open Summary file and check the peak list files, upl, and ovw to determine the number of assignments and violations and write out the summary file 

tpeak,tsingle,tamb,tnotused,tnota,tdia,tincr,tupl,tviol  = 0, 0, 0, 0, 0, 0, 0, 0, 0

for x in range(len(cya_plists)):
  plistn = cya_plists[x]
  plist = cya_plists[x]
  upl_re_str = r'.*plist ' + str(x+1) + '([^0-9]|$)'
  upl = [line.strip() for line in open(fupl).readlines() if line.strip() and re.match(upl_re_str, line)]
  tupl = tupl + len(upl)
  viol_re_str = r'.*list ' + str(x+1) + '([^0-9]|$)'
  viol = [line.strip() for line in open(fovw).readlines() if line.strip() and re.match(viol_re_str, line)]
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
  peaklisttxt = peaklisttxt + "{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%\n".format(plist,linepad,peak,len(upl),len(viol), single,amb,notused,nota,dia,incr, 100*((single+amb+dia)/peak))
  shortsumtext = shortsumtext + "{:<}{:} {:^6d}   {:^3.1f}% ({:})\n".format(plist,linepad,peak,100*((single+amb+dia)/peak),(single+amb+dia))
  # print("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%".format(plist,linepad, peak,len(upl),len(viol),single,amb,notused,nota,dia,incr, 100*((single+amb+dia)/peak)))
peaklisttxt = peaklisttxt + "{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%\n".format('Total',pad[5:],tpeak,tupl,tviol,tsingle,tamb,tnotused,tnota,tdia,tincr,100*((tsingle+ tamb+ tdia)/tpeak))
shortsumtext = shortsumtext + "{:<}{:} {:^6d}   {:^3.1f}% ({:})\n\n".format('Total',pad[5:],tpeak,100*((tsingle+ tamb+ tdia)/tpeak),(tsingle+ tamb+ tdia))
#print("{:<}{:} {:^6d}  {:^4d} {:^4d} {:^7d} {:^9d} {:^7d} {:^5d}  {:^8d} {:^13}  {:^3.1f}%".format('Total',pad[5:],tpeak,tupl,tviol,tsingle,tamb,tnotused,tnota,tdia,tincr, 100*((tsingle+ tamb+ tdia)/tpeak)))
print(peaklisttxt)
checkcons.write(peaklisttxt +'\n\n')
shortsum.write(shortsumtext)

outpml = open(outdir + pdbname + '.pml','w')
outpml.write('load ./' + in_pdb+'\n')
outpml.write('set dash_gap, 0.05\n')
outpml.write('set_color royalblue = [65,105,225]\nset_color forest = [34,139,34]\nset_color yellowgreen = [154,205,50]\nset_color darkorange = [255,140,0]\nset_color purple = [128,0,128]\nset_color lightseagreen = [32,178,170]\nset_color darkkhaki = [189,183,107]\nset_color peru = [205,133,63]\nset_color saddlebrown = [139,69,19]\nset_color gold = [255,215,0]\nset_color navy = [0,0,128]\nset_color darkturquoise = [0,206,209]\nset_color pink = [255,192,203]\nset_color cyan = [0,255,255]\nset_color paleturquoise = [175,238,238]\nset_color lightsalmon = [255,160,122]\nset_color khaki = [240,230,140]\nset_color yellowgreen = [154,205,50]\nset_color thistle = [216,191,216]\nset_color aquamarine = [127,255,212]\nset_color plum = [221,160,221]\nset_color lightpink = [255,182,193]\nset_color mediumvioletred = [199,21,133]\nset_color firebrick = [178,34,34]\nset_color lightcoral = [240,128,128]\nset_color deeppink = [255,20,147]\nset_color hotpink = [255,105,180]\nset_color mediumpurple = [147,112,219]\nset_color navy = [0,0,128]\nset_color cornflowerblue = [100,149,237]\n')
outpml.write('color gray60, all\n')
outpml.write('show sticks, {:} and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\n hide sticks, elem H\nhide sticks, name N+C\n'.format(pdbname))
outpml.write('color paleturquoise, {:} and resn ILE\ncolor lightsalmon, {:} and resn LEU\ncolor khaki, {:} and resn VAL\ncolor yellowgreen, {:} and resn ALA\ncolor thistle, {:} and resn MET\ncolor aquamarine, {:} and resn THR\ncolor lightpink, {:} and resn TYR\ncolor plum, {:} and resn PHE\n'.format(pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname))
outpml.write('color gold, elem S\ncolor red, elem O\ncolor blue, elem N\n')
outpml.write('split_states ' + pdbname + '\n')
for y in range(2,21,1):
  outpml.write('align {:}_{:04d}, {:}_0001\n'.format(pdbname,y, pdbname))
outcmx = open(outdir + pdbname + '.cxc','w')
outcmx.write('open ../'+ in_pdb+'\n')
outcmx.write('color #1 gray(150)\n')
outcmx.write('match #1.2-20 to #1.1\n')
outcmx.write('color #1:ile paleturquoise target a\ncolor #1:leu lightsalmon  target a\ncolor #1:val khaki target a\ncolor #1:ala yellowgreen  target a\ncolor #1:met thistle target a\ncolor #1:thr aquamarine target a\ncolor #1:phe plum target a\ncolor #1:tyr lightpink target a\n')
outcmx.write('show #1:thr,met,ala,leu,val,ile,phe,tyr\ncartoon suppress false\n')
outcmx.write('label #1.1 text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\n''label ontop false\n')
outcmx.write('ui tool show "Side View"\n#ui mousemode right distance\n')

## ---------------------------------------------------------------------------
## Creating the pseudobond files and group strings for rendering the constraints in chimera/pymol
## Make pseudo bond and gorup statements for poor, long, and short upl entries from designated upl file 
for conect in ConectionTypes:
  MasterDict["{:}_pb".format(conect)] = []
  MasterDict["group{:}".format(conect)] = "'group {:}, '".format(conect)
poorpbout = open(outdir +'pseudobonds/' + outname + '_poor.pb','w')
poorpbout.write("; halfbond = false\n; color = mediumvioletred\n; radius = 0.1\n; dashes = 0\n")
longpbout = open(outdir +'pseudobonds/' + outname + '_long.pb','w')
longpbout.write("; halfbond = false\n; color = firebrick\n; radius = 0.1\n; dashes = 0\n")
pdiffpbout = open(outdir +'pseudobonds/' + outname + '_pdiff.pb','w')
pdiffpbout.write("; halfbond = false\n; color = firebrick\n; radius = 0.1\n; dashes = 0\n")
shortpbout = open(outdir +'pseudobonds/' + outname + '_short.pb','w')
shortpbout.write("; halfbond = false\n; color = light coral\n; radius = 0.1\n; dashes = 0\n")
uviolpbout = open(outdir +'pseudobonds/' + outname + '_viol_upls.pb','w')
uviolpbout.write("; halfbond = false\n; color = hotpink\n; radius = 0.1\n; dashes = 0\n")
finalupls = [["###Violated Restraints\n"],["###Poor/Low Support\n"],["###Long Distance Restraints (d >= 6.0)\n"],["###Short Distance Restraints (d <= 3.0)\n"],["###Good Restraints\n"]]

# cmxphisel, cmxchisel, cmxphiviol, cmxchiviol = 'name phipsisel #angmn:', 'name chisel #angmn:', 'name phipsiviol #angmn:', 'name chiviol #angmn:'
# pmlphisel, pmlchisel, pmlphiviol, pmlchiviol = 'phipsi and resi ','chi and resi ', 'phipsi and resi ', 'chi and resi '

### Go through the final overview file and extract information about violated distance and angle restraints 
Filtered = []
violdict, Upperdict, Lowerdict, dihedviol = {}, {}, {}, {}
viol_upls= 'group viol_upl, ';MasterDict['viol_upls'] = viol_upls
phiviol,chiviol,violpeaks,diheviols =  [], [], [], []
v = 0
vphicount, vpsicount,vchi1count,vchi2count, vothercount = 0,0,0,0,0
hbondline = ''
for line in open(fovw).readlines():
  if line[4:9] == 'Upper' or line[4:9] == 'Lower':
    dviol = line.split()
    if dviol[2]+dviol[1] not in Ambiguous.keys() and dviol[6]+dviol[5] not in Ambiguous.keys(): ## skip the ambiguous entries they are not used in refinement 
      atom1 = dviol[1]
      if dviol[2]+dviol[1] in replacements.keys():
        atom1 = replacements[dviol[2]+dviol[1]]
      atom2 = dviol[5]
      if dviol[6]+dviol[5] in replacements.keys():
        atom2 = replacements[dviol[6]+dviol[5]]
      cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[2]],dviol[3],dviol[1],AAA_dict[dviol[6]],dviol[7],dviol[5])
      cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[6]],dviol[7],dviol[5],AAA_dict[dviol[2]],dviol[3],dviol[1])
      if int(dviol[9]) >= 6:
        pbout = uviolpbout
        grpstr = "uplviol"
        outline = ' #viol in {:} by +{:}\n'.format(dviol[9],dviol[10])
        if 'peak' not in line:
          if line[4:9] == 'Upper':
            if 'O' not in dviol[1]:
              upldf.loc[AAA_dict[dviol[6]] + dviol[7],'viol input'] = upldf.loc[AAA_dict[dviol[6]] + dviol[7],'viol input'] + 1
              upldf.loc[AAA_dict[dviol[2]] + dviol[3],'viol input'] = upldf.loc[AAA_dict[dviol[2]] + dviol[3],'viol input'] + 1
              Upperdict[cons1] = outline
              Upperdict[cons2] = outline
          if line[4:9] == 'Lower':
            outline = ' #viol in {:} by -{:}\n'.format(dviol[9],dviol[10])
            Lowerdict[cons1] = outline
            Lowerdict[cons2] = outline
        if 'peak' in line and 'list' in line:
          violdict[cons1] = outline
          violdict[cons2] = outline
          for line2 in open(fupl).readlines():
            cns = line2.split()
            if cns[8] == line[90:].split()[1] and cns[10] == line[90:].split()[3] and cns[2] == dviol[1] and cns[5] == dviol[5]:
              if cns[0] != cns[3] and '#SUP' in line2:
                upldf.loc[AAA_dict[cns[1]] + cns[0],'viol'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'viol'] + 1
                upldf.loc[AAA_dict[cns[4]] + cns[3],'viol'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'viol'] + 1
                violpeaks.append(cons1)
                finalupls[0].append(line2)
                Filtered.append(line2)
        v+=1
        pbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(dviol[3], atom1, dviol[7],atom2))
        outpml.write('distance viol{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(v), pdbname, dviol[3], atom1, pdbname, dviol[7], atom2))
        viol_upls = viol_upls + "viol"+str(v) + ' '
  if line[4:9] == 'Angle':
    dang = line.split()
    dihedviol[AAA_dict[dang[2]] + dang[3] + dang[1].replace('CHI21','CHI2')] = r'$\{:}$ viol in {:} by {:}'.format(dang[1].lower(), dang[6], dang[7])
    diheviols.append(r'$\{:}$ viol in {:} by {:}'.format(dang[1].lower(), dang[6], dang[7]))
    angle = dang[1].replace('CHI21','CHI2')
    try:
      exec('v{:}count = v{:}count + 1'.format(angle.lower(),angle.lower()))
    except NameError:
      vothercount+=1
    if dang[1] == 'PHI' or dang[1] == 'PSI':
      if dang[3] not in phiviol:
        phiviol.append(dang[3])
        # cmxphiviol = cmxphiviol + dang[3] + ','
        # pmlphiviol = pmlphiviol + dang[3] + '+'
    if 'CHI' in dang[1]:
      if dang[3] not in chiviol:
        chiviol.append(dang[3])
        # cmxchiviol = cmxchiviol + dang[3] + ','
        # pmlchiviol = pmlchiviol + dang[3] + '+'
  if 'Hydrogen bonds in 6' in line:
    hbondline = line

usedupls,qupldict, upldict, upldict2 = {}, {}, {}, []
finalupl,poorcons2, show, shortcons2,longcons2,diffcons2,sidelist = [],[],[],[],[],[],[]
poorcons, shortcons, longcons, pdiffcons = 'group poor, ', 'group short, ', 'group long, ', 'group pdiffusion, '
MasterDict['poorcons'] = poorcons; MasterDict['shortcons'] = shortcons;MasterDict['longcons']=longcons; MasterDict['pdiffcons'] =pdiffcons
for line in open(fupl).readlines():
  if line.split() and '#SUP' in line:
    cns = line.split()
    cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
    cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5],AAA_dict[cns[1]],cns[0],cns[2])
    if cons1 not in upldict.keys():
      upldict[cons1] = "{:3.2f}".format(float(cns[6]))
    if cons2 not in upldict.keys():
      upldict[cons2] = "{:3.2f}".format(float(cns[6]))
    upldict2.append("{:} peak {:} from {:}".format(cons1,cns[8],cya_plists[int(cns[10])-1]))
    upldict2.append("{:} peak {:} from {:}".format(cons2,cns[8],cya_plists[int(cns[10])-1]))
fovwlines = open(fovw).readlines().index(hbondline) + 2 

for line in open(fovw).readlines()[fovwlines:]:
  if line.strip():
    if line.strip()[0] == 'H':
      hbc = line.split()
      if len(hbc[0]) <= 3: doner = hbc[0].replace('H','N')
      if len(hbc[0]) == 4: doner = hbc[0][:-1].replace('H','N')
      usedupls['{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[hbc[1]],hbc[2],hbc[0],AAA_dict[hbc[5]],hbc[6],hbc[4])] = ' # found in {:} structures\n'.format(hbc[7])
      usedupls['{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[hbc[5]],hbc[6],hbc[4],AAA_dict[hbc[1]],hbc[2],hbc[0])] = ' # found in {:} structures\n'.format(hbc[7])
      usedupls['{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[hbc[1]],hbc[2],doner,AAA_dict[hbc[5]],hbc[6],hbc[4])] = ' # found in {:} structures\n'.format(hbc[7])
      usedupls['{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[hbc[5]],hbc[6],hbc[4],AAA_dict[hbc[1]],hbc[2],doner)] = ' # found in {:} structures\n'.format(hbc[7])

ADpairs,Assignments = [],[]
noalines = open(noa).readlines()

Swapped = {}
log = glob.glob(os.path.join(cwd + 'log*'))[0]
for line in open(log).readlines():
  if line.strip() and line.strip().split()[-1] == 'swapped':
    res = line.strip().split()
    group1 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[2])
    group2 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[3])
    Swapped[group1] = group2
    Swapped[group2] = group1
for x in range(len(noalines)):
  if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
    if '0 out of 0' not in noalines[x+1] and 'diagonal' not in noalines[x]:
      index = 0
      if '0 out of' in noalines[x+1]: 
        index = 3
      for y in range(2,int(noalines[x+1].split()[index])+2,1):
        cns = noalines[x+y].strip().split()
        if not re.match(r'^\s+[A-Z]+',noalines[x+y]):
          atom1,resn1, resi1, atom2, resn2, resi2 = cns[1],cns[2],int(cns[3]), cns[5], cns[6], int(cns[7])
        if re.match(r'^\s+[A-Z]+',noalines[x+y]):
          atom1,resn1, resi1, atom2, resn2, resi2 = cns[0],cns[1],int(cns[2]), cns[4], cns[5], int(cns[6])
        group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
        group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
        if group1 in Swapped.keys(): group1 = Swapped[group1]
        if group2 in Swapped.keys(): group2 = Swapped[group2]
        if '{:}-{:}'.format(group1,group2) not in ADpairs or '{:}-{:}'.format(group2,group1) not in ADpairs:
          ADpairs.append('{:}-{:}'.format(group1,group2))
        if group1 not in Assignments: Assignments.append(group1)
        if group2 not in Assignments: Assignments.append(group2)

Assignments = sorted(Assignments , key = lambda x:(x.split('-')[0][1:], x.split('-')[1]))
posdiffout = open(outdir + 'probable_diffusion_log.txt','w')
import GetDistances as distance
distDF,fillteredDF = distance.examin(in_pdb,ADpairs,Assignments,fupl)
# exit()
i = 1
im = 0
for line in open(fupl).readlines():
  if line not in Filtered:
    if '#SUP' not in line: ## exclude ambiguous restraints
      pass 
    else:
      cns = line.split()
      if cns[0] == cns[3]: ## exclude intramolecular restraints
        pass
      if cns[0] != cns[3]:
        if cns[1]+cns[2] in ConTypeDict.keys():ct1 = ConTypeDict[cns[1]+cns[2]]
        if cns[1]+cns[2] not in ConTypeDict.keys(): ct1 = 'Other'
        if cns[4]+cns[5] in ConTypeDict.keys():ct2 = ConTypeDict[cns[4]+cns[5]]
        if cns[4]+cns[5] not in ConTypeDict.keys(): ct2 = 'Other'
        ctype = ConTypeDict["{:}-{:}".format(ct1,ct2)]
        pblist = MasterDict[ctype + '_pb']
        atom1 = cns[2]
        atom2 = cns[5]
        if atom1 == 'H': atom1 = 'N'
        if atom2 == 'H': atom2 = 'N'
        if cns[1]+cns[2] in replacements.keys():
          atom1 = atom1.replace(cns[2], replacements[cns[1]+cns[2]])
        if cns[4]+cns[5] in replacements.keys():
          atom2 = atom2.replace(cns[5], replacements[cns[4]+cns[5]])
        upldf.loc[AAA_dict[cns[1]] + cns[0],'cya'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'cya'] + 1
        upldf.loc[AAA_dict[cns[4]] + cns[3],'cya'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'cya'] + 1
        group1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
        group2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
        outline = ' #{:3.2f} #peak {:} #plist {:}\n'.format(float(cns[6]),cns[8],cns[10])
        cons1 = '{:}-{:}'.format(group1,group2)
        usedupls['{:}-{:}'.format(group1,group2)] = outline
        usedupls['{:}-{:}'.format(group2,group1)] = outline
        ## Traslate amides H to N to search input upls but don't mess up the connections drawn in pymol/chimera
        atoms1 = atom1
        atoms2 = atom2
        if atom1 == 'H': atoms1 = 'N'
        if atom2 == 'H': atoms2 = 'N'
        cons1a = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],atoms1,AAA_dict[cns[4]],cns[3],atoms2)
        cons2a = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],atoms2,AAA_dict[cns[1]],cns[0],atoms1)
        usedupls[cons1a] = outline
        usedupls[cons2a] = outline
        if (float(cns[0]) > 1000 and float(cns[3]) < 1000) or (float(cns[0]) < 1000 and float(cns[3]) > 1000):
          im+=1
          ctype = 'Interamolecular'
          pblist = intermolecular_pb
          outpml.write('distance intermol{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(im), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
          exec('groupintermolecular = groupintermolecular+ "intermo{:} "'.format(im))
          pblist.append('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
        ## Make sure all connections to side chains are shown
        if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
          sidelist.append(cns[0])
        if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
          sidelist.append(cns[3])
        try:
          d = float(distDF.loc[group1,group2].split()[0])
        except (AttributeError, KeyError):
          d = 0
        if d >= 7.7 and line not in Filtered:
          common = fillteredDF.dropna(subset=[group1,group2]).index.tolist()
          if len(common) >= 2:
            common.extend([group1,group2])
            upldf.loc[AAA_dict[cns[1]] + cns[0],'long'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'long'] + 1
            upldf.loc[AAA_dict[cns[4]] + cns[3],'long'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'long'] + 1
            qupldict['{:}-{:}'.format(group1,group2)] = "prob diff {:}A ".format(distDF.loc[group1,group2])
            qupldict['{:}-{:}'.format(group2,group1)] = "prob diff {:}A ".format(distDF.loc[group1,group2])
            i+=1
            pdiffpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
            diffcons2.append('{:}-{:}'.format(group1,group2))
            outpml.write('distance pdiffusion{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
            pdiffcons = pdiffcons + 'pdiffusion{:} '.format(i)
            Filtered.append(line)
            posdiffout.write('{:}-{:}\n'.format(group1,group2))
            posdiffout.write('long distance upl {:} heavy {:}\n'.format(cns[6],d))
            posdiffout.write(fillteredDF.loc[common,common].to_string())
            posdiffout.write('\n')
            Filtered.append(line)
        if line not in Filtered and float(cns[12]) < 0.5:
          i+=1
          poorpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
          poorcons2.append('{:}-{:}'.format(group1,group2))
          outpml.write('distance poor{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
          poorcons = poorcons + 'poor{:} '.format(i)
          Filtered.append(line)
          finalupls[1].append(line)
        if line not in Filtered and float(cns[12]) > 0.5:
          if float(cns[6]) >= 6.0:
            upldf.loc[AAA_dict[cns[1]] + cns[0],'long'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'long'] + 1
            upldf.loc[AAA_dict[cns[4]] + cns[3],'long'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'long'] + 1
            qupldict['{:}-{:}'.format(group1,group2)] = "long "
            qupldict['{:}-{:}'.format(group1,group2)] = "long "
            i+=1
            longpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
            longcons2.append('{:}-{:}'.format(group1,group2))
            outpml.write('distance long{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
            longcons = longcons + 'long{:} '.format(i)
            Filtered.append(line)
            finalupls[2].append(line)
          if float(cns[6]) <= 3.00:
            if abs(int(cns[0])- int(cns[3])) > 1:
              i+=1
              shortpbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
              shortcons2.append('{:}-{:}'.format(group1,group2))
              outpml.write('distance short{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
              qupldict['{:}-{:}'.format(group1,group2)] = "short "
              qupldict['{:}-{:}'.format(group1,group2)] = "short "
              shortcons = shortcons + 'short{:} '.format(i)
              finalupls[3].append(line)
              Filtered.append(line)
          if float(cns[6]) > 3.00 and float(cns[6]) < 6.0:
            i+=1
            outpml.write('distance UPL{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
            pblist.append('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
            MasterDict['group' + ctype] = MasterDict['group' + ctype] + '+ "UPL{:} "'.format(i)
            Filtered.append(line)
            finalupls[4].append(line)
if im != 0: 
  print('found {:} intermolecuar connections'.format(im))
poorpbout.close()
longpbout.close()
shortpbout.close()
uviolpbout.close()

# ---------------------------------------------------------------------------
# Run the cycle7.noa analysis and generate peak list identifying peaks as 
# unused, no assignment, and questionable

assigndict = noaa.analize_noa(Seqdict, violdict, qupldict, upldict, pad, upldict2, distDF, fillteredDF ,posdiffout)

mn = 1
for x in range(len(ConectionTypes)):
  pbs = MasterDict['{:}_pb'.format(ConectionTypes[x])]
  if len(pbs) > 1:
    mn+=1
    pbout = open('{:}pseudobonds/{:}_{:}.pb'.format(outdir, outname, ConectionTypes[x]),'w')
    pbout.write("; halfbond = false\n; color = " + colors[x] + "\n; radius = 0.1\n; dashes = 0\n")
    pbout.writelines(pbs)
    outcmx.write('open pseudobonds/{:}_{:}.pb\n'.format(outname, ConectionTypes[x]))
    # outcmx.write('color #{:} {:}\n'.format(str(mn),colors[x]))
    groupstr = MasterDict['group' + ConectionTypes[x]]
    outpml.write(groupstr + '\n')
    outpml.write('color {:}, {:}\n'.format(colors[x],ConectionTypes[x]))

for (group, color) in [('poor','mediumvioletred'),('long','firebrick'), ('pdiff','firebrick'), ('short', 'lightcoral'),('viol_upls', 'deeppink')]:
  mn+=1
  outcmx.write('open pseudobonds/' + outname + '_' + group + '.pb\n')
  # outcmx.write('color #{:} {:}\n'.format(str(mn),color))
for (group,gname, color) in [('poorcons','poor','mediumvioletred'),('longcons','long','firebrick'), ('pdiffcons','pdiffusion','firebrick'),('shortcons','short', 'lightcoral'),('viol_upls','viol_upl', 'deeppink')]:
  grpstr = MasterDict[group]
  outpml.write(grpstr + '\n')
  outpml.write('color {:}, {:}\n'.format(color, gname))

#### Write out the filtered upl list, which does not contain ambiguous (QQ) restraints 
#### and has sorted the restraints into 5 labeled categories
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
  outpb.write("; halfbond = false\n; color = cyan\n; radius = 0.15\n; dashes = 10\n")
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
      # print(line.rstrip())
      upldf.loc[AAA_dict[cns[1]] + cns[0],'input'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'input'] + 1
      upldf.loc[AAA_dict[cns[4]] + cns[3],'input'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'input'] + 1
      cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
      for atom1 in atoms1:
        for atom2 in atoms2:
          u+=1
          outpml.write('distance {:}{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(uplfile.replace('.upl',''),str(u), pdbname, cns[0], atom1, pdbname, cns[3], atom2))
          pmlgroup = pmlgroup + uplfile.replace('.upl','') + str(u) + ' '
          if 'missing' in line:
            outpb.write('#1.1:{:}@{:} #1.1:{:}@{:} blue\n'.format(cns[0], atom1, cns[3],atom2))
          if cons in Upperdict.keys():
            outpb.write('#1.1:{:}@{:} #1.1:{:}@{:} hotpink\n'.format(cns[0], atom1, cns[3],atom2))
          else:
            outpb.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], atom1, cns[3],atom2))
      if (cns[1] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[2] not in ['N','H']) and cns[0] not in sidelist:
        sidelist.append(cns[0])
      if (cns[4] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[3] not in sidelist:
        sidelist.append(cns[3])
  outpml.write(pmlgroup + '\n')
  outpml.write('color cyan,' + uplfile.replace('.upl','') + '\n')
  mn+=1
  outcmx.write('open pseudobonds/' + uplfile.replace('.upl','.pb') + '\n')
selhbond = 'name hbond #1.1:'
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
      if (cns[0],cns[3],cns[5]) not in hbonsl:
        h+=1 
        hbonsl.append((cns[0],cns[3],cns[5].replace('N','H')))
        cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
        if cons in Upperdict.keys():
          hbond.write('#1.1:{:}@{:} #1.1:{:}@{:} hotpink\n'.format(cns[0], cns[2], cns[3],cns[5]))
        else:
          hbond.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[0], cns[2], cns[3],cns[5]))
        outpml.write('distance hbond{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(h), pdbname, cns[0], cns[2].replace('H','N'), pdbname, cns[3], cns[5].replace('H','N')))
        hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
        if cns[0] not in selhbond:
          selhbond = selhbond +'{:},'.format(cns[0])
        if cns[3] not in selhbond:
          selhbond = selhbond +'{:},'.format(cns[3])
        if cns[2] not in ['O','H','N']: sidelist.append(cns[0])
        if cns[5] not in ['O','H','N']: sidelist.append(cns[3])
hbond.close()
sidechains = 'name sidechains #1:'
for res in sidelist:
  sidechains = sidechains + res + ','
outcmx.write(sidechains[:-1] + '\nshow sidechains\n')
outpml.write('show sticks, resi '+ sidechains.split(':')[-1][:-1].replace(',','+')+'\n')
outcmx.write('hide H\n''show #1.1@H,N,C target a\n')
outpml.write(hbgroupline + '\n')
outpml.write('color pink, hbond\n')
outpml.write('show sticks,name N+O+C+H and resi '+ selhbond.split(':')[-1][:-1].replace(',','+')+'\n')
selhbond = selhbond[:-1] + '@O,N,C\nshow hbond target a\n'
outcmx.write('open pseudobonds/' + 'hbond.pb\n')
outcmx.write(selhbond)
upls.extend(hbonds)
mn+=1

## Examine input upl files and update lines for entries which are violated 10 or more times, and identify proton-proton restraints that support heavy atoms based restraints 

for uplfile in upls:
  found_upls, tupls, mupls, vupls  = 0, 0, 0, 0
  newlines = []
  violupl = open(outdir + uplfile.replace('.upl','_viol.upl'),'w')
  matchedupl = open(outdir + uplfile.replace('.upl','_found.upl'),'w')
  for line in open(uplfile).readlines():
    newline = ''
    if line.strip():
      if '#' not in line.split()[0]:
        if 'missing' not in line: tupls+=1
        if 'missing' in line: mupls+=1
        cns = line.split()
        atoms1 = [cns[2]]
        atoms2 = [cns[5]]
        cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
        if cns[1]+cns[2] in Ambiguous.keys(): atoms1 = Ambiguous[cns[1]+cns[2]].split(',')
        if cns[4]+cns[5] in Ambiguous.keys(): atoms2 = Ambiguous[cns[4]+cns[5]].split(',')
        for atom1 in atoms1:
          for atom2 in atoms2: 
            cons = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],atom1,AAA_dict[cns[4]],cns[3],atom2)
            if cons in Upperdict.keys():
              vupls+= 1
              newline = line.replace('\n', Upperdict[cons])
              violupl.write(newline)
            if cons in usedupls.keys():
              found_upls+= 1
              if 'hbond' not in uplfile:
                upldf.loc[AAA_dict[cns[1]] + cns[0],'found input'] = upldf.loc[AAA_dict[cns[1]] + cns[0],'found input'] + 1
                upldf.loc[AAA_dict[cns[4]] + cns[3],'found input'] = upldf.loc[AAA_dict[cns[4]] + cns[3],'found input'] + 1
              if len(newline) > 1: newline = newline.replace('\n', usedupls[cons])
              if len(newline) < 1: newline = line.replace('\n', usedupls[cons])
              matchedline = line.replace('\n', usedupls[cons])
              matchedupl.write(matchedline)
    if len(newline) < 1:
      newline = line
    newlines.append(newline)
  checkcons.write('{:} input upls from {:}\n'.format(mupls + tupls, uplfile))
  shortsum.write('{:} input upls from {:}\n'.format(mupls + tupls, uplfile))
  checkcons.write('    {:} of {:} of assignable input upls found\n'.format(found_upls,tupls))
  # shortsum.write('    {:} of {:} of assignable input upls found\n'.format(found_upls,tupls))
  checkcons.write('    {:} of upls violated in 10 or more structures\n'.format(vupls))
  shortsum.write('    {:} of upls violated in 10 or more structures\n'.format(vupls))
  if mupls >=1:
    checkcons.write('    {:} of upls missing assignment\n'.format(mupls))
    # shortsum.write('    {:} of upls missing assignment\n'.format(mupls))

  fout = open(outdir + uplfile,'w')
  fout.writelines(newlines)
  fout.close()
  violupl.close()
  matchedupl.close()
shortsum.write('\n')
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
checkcons.write('\n\n')
phiaco, chiaco = [], []
phicount, psicount,chi1count,chi2count, othercount, total = 0,0,0,0,0,0
phipsidict,chidict,plotdict = {}, {}, {}
for aco in dihed:
  for line in open(aco):
    if line.split():
      if re.match(r'^\s*#', line):
        continue
      else:
        total+=1
        ang = line.split()
        angle = ang[2].replace('CHI21','CHI2')
        try:
          exec('{:}count = {:}count + 1'.format(angle.lower(),angle.lower()))
        except NameError:
          othercount+=1
        ang_min = float(ang[3])
        ang_max = float(ang[4])
        if 'P' in ang[2]:
          angdict = phipsidict
          if ang[0] not in phiaco:
            phiaco.append(ang[0])
            # cmxphisel = cmxphisel + ang[0] + ','
            # pmlphisel = pmlphisel + ang[0] + '+'
        if 'CHI' in ang[2]:
          angdict = chidict
          if ang[0] not in chiaco:
            chiaco.append(ang[0])
          if ang_min < 0:
            ang_min = ang_min + 360.0
            ang_max = ang_max + 360.0
            # cmxchisel = cmxchisel + ang[0] + ','
            # pmlchisel = pmlchisel  + ang[0] + '+'
        outline = r"$\{:}$  {:} - {:}".format(angle.lower(), ang_min,ang_max)
        plotdict[AAA_dict[ang[1]] + ang[0] + angle] = [ang_min, ang_max]
        if AAA_dict[ang[1]] + ang[0] not in angdict.keys():
          angdict[AAA_dict[ang[1]] + ang[0]] = [[outline,'black']]
        else: 
          angdict[AAA_dict[ang[1]] + ang[0]].append([outline,'black'])

angle_text = "Total of {:} dihedral restraints:\n       input viol\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n\n".format(total, 'Phi', phicount, vphicount, 'Psi', psicount, vpsicount , 'Chi1', chi1count ,vchi1count, 'Chi2', chi2count, vchi2count)

## ---------------------------------------------------------------------------
## Run the GetDihed.py to determine phi, psi, chi1 and chi2 and plot them
## for all 20 structures
## ---------------------------------------------------------------------------
print('Extracting dihedrals')
text = [['CYANA UPL','#9acd32'],['long UPL','#800080'],['Violated UPL','#ffa500'],['Input UPL','#6495ed'],['Found Input UPL','navy'],['Violate Input UPL','#db7093']]
DAramalist, DArotalist = Dihed.extract(in_pdb, ASequence, outdir, upldf, phipsidict, chidict, plotdict,dihedviol,text)
print('finished plotting dihedrals')
armn = mn
routmn = mn+1
mn+=1
outcmx.write('open ../{:} maxModels 1\nrename #{:} angle_restraints\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(in_pdb,armn,armn,armn))
outcmx.write('open ../{:} maxModels 1\nrename #{:} rama_outliers\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(in_pdb,routmn,routmn,routmn))
ramalist, rotalist = [],[]
for line in DAramalist:ramalist.append(line.split()[0][1:])
for line in DArotalist:rotalist.append(line.split()[0][1:])
anglesout = [[phiaco, "phipsi","purple",armn], [chiaco,"chi1chi2","navy",armn], 
[phiviol,"viol_phipsi","mediumpurple",armn], [chiviol,"chi1chi2","cornflowerblue",armn],
[ramalist,"dissallowed_phipsi","mediumvioletred",routmn], [rotalist,"dissallowed_chi1chi2","medium violet red",routmn]]
for listn, aname, color, modle in anglesout:
  outpml.write('create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\n'.format(aname,pdbname,aname,aname))
  rline = ""
  if len(listn) > 1:
    for resn in listn:
      rline = rline + "{:},".format(resn)
    rline = rline[:-1]
    if 'phi' in aname:
      outcmx.write("name {:} #{:}:{:}\ncolor {:} {:} target c\n".format(aname,modle,rline,aname,color))
      outpml.write("create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\ncolor {:}, {:} and resn {:}\n".format(aname,pdbname,aname,aname,color,aname,rline.replace(',','+')))
    if 'chi' in aname:
      outcmx.write("name {:} #{:}:{:}\nshow {:} target a\ncolor {:} {:} target a\nhide #{:}@H*,N,O target a\ncolor byhetero target a\n".format(aname,modle,rline,aname,aname,color,modle))
      outpml.write("create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\nshow sticks, {:} and resn {:}\ncolor {:}, {:} and resn {:}\n".format(aname,pdbname,aname,aname,aname,rline.replace(',','+'),color,aname,rline.replace(',','+')))

print(angle_text[:-2])
checkcons.write(angle_text)
shortsum.write(angle_text)
overviewtxt = '{:3.0f} Residues with Disallowed Phi/Psi \n{:3.0f} Residues with Disallowed Chi1/Chi2\n{:3.0f} Violated Distance Restraints\n{:3.0f} Low Support Restraints\n{:3.0f} Probable Diffusion Distance Restraints\n{:3.0f} Long Distance Restraints d >= 6.0\n{:3.0f} Short Distance Restraints d <= 3.0\n\n'.format(len(DAramalist),len(DArotalist),len(violpeaks),len(poorcons2),len(diffcons2),len(longcons2),len(shortcons2))
print(overviewtxt)
checkcons.write(overviewtxt)
bad_entries = []
checkcons.write('### {:3.0f}  Violated Distance Restraints ###\n'.format(len(violpeaks)))
# violpeaks = sorted(violpeaks, key = lambda x: (x.split()[10],x.split()[8]))
violpeaks = sorted(violpeaks, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
for viol in violpeaks:
  if viol in assigndict.keys():
    try:
      checkcons.write('{:}  {:3.2f}A ({:}): {:}'.format(viol,float(upldict[viol]),len(assigndict[viol]), violdict[viol]))
    except KeyError:
      print('Found Error {:} violation not in upl list'.format(viol))
    # checkcons.write('{:}  {:3.2f}A ({:}): {:}'.format(viol,float(upldict[viol]),len(assigndict[viol]), violdict[viol]))
    checkcons.writelines(assigndict[viol])
    checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Poor/Low Support constraints to the summary file
poorcons2 = sorted(poorcons2, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
checkcons.write('### {:3.0f} Low Support Restraints ###\n'.format(len(poorcons2)))
for con in poorcons2:
  if con in assigndict.keys():
    checkcons.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
    checkcons.writelines(assigndict[con])
    checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Probable Diffusion Distance constraints to the summary file
diffcons2 = sorted(diffcons2, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
checkcons.write('### {:3.0f} Probable Diffusion Distance Restraints ###\n'.format(len(diffcons2)))
for con in diffcons2:
  if con in assigndict.keys():
    checkcons.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
    checkcons.writelines(assigndict[con])
    checkcons.write('\n')
#### Write out Long Distance constraints to the summary file
longcons2 = sorted(longcons2, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
checkcons.write('### {:3.0f} Long Distance Restraints d >= 6.0\u00C5 ###\n'.format(len(longcons2)))
for con in longcons2:
  if con in assigndict.keys():
    checkcons.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
    checkcons.writelines(assigndict[con])
    checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Short Distance constraints to the summary file
checkcons.write('### {:3.0f} Short Distance Restraints d <= 3.0\u00C5 ###\n'.format(len(shortcons2)))
shortcons2 = sorted(shortcons2, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
for con in shortcons2:
  if con in assigndict.keys():
    checkcons.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
    checkcons.writelines(assigndict[con])
    checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Violated Dihedral Restraints to the summary file
checkcons.write('### {:3.0f} Violated Dihedral Restraints ###\n'.format(len(diheviols)))
diheviols = sorted(diheviols, key = lambda x: (x.split()[0][1:], x.split()[1]))
for con in diheviols:
  checkcons.writelines(con)
checkcons.write('\n\n')
#### Write out Residues with Disallowed Phi/Psi to the summary file
checkcons.write('### {:3.0f} Disallowed Phi/Psi Dihedral Restraints ###\n'.format(len(DAramalist)))
for con in DAramalist:
  checkcons.writelines(con)
checkcons.write('\n\n')
#### Write out Residues with Disallowed Chi1/Chi2 to the summary file
checkcons.write('### {:3.0f} Disallowed Chi1/Chi2 Dihedral ###\n'.format(len(DArotalist)))
for con in DArotalist:
  checkcons.writelines(con)
checkcons.write('\n\n')
checkcons.close()

for line in open(fovw).readlines():
  if line.strip():
    if line.split()[0] == 'Ave':
      shortsum.write("Target Function {:}\n".format(line.split()[1]))
    if line.split()[0] == 'Average':
      shortsum.write(line.strip()[8:52] + '\n')
shortsum.close()

## ---------------------------------------------------------------------------
### Creating Model coloring residues based on the number of NOE restraints 
mn+=1
outcmx.write('open ../' + in_pdb+ ' maxModels 1\nrename #{:} noes\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(mn,mn,mn))
outcmx.write('color name c0 rgb(255,205,0)\ncolor name c2 rgb(156,217,59)\ncolor name c4 rgb(52,182,121)\ncolor name c6 rgb(42,117,142)\ncolor name c8 rgb(59,81,139)\ncolor name c10 rgb(20,64,110)\n')
outpml.write('set_color c0 = [255,205,0]\nset_color c2  = [156,217,59]\nset_color c4  = [52,182,121]\nset_color c6  = [42,117,142]\nset_color c8  = [59,81,139]\nset_color c10  = [20,64,110]\n')
outpml.write('create noes, {:}_0001\ncolor gray60,phi-psi\nhide sticks, noes\n'.format(pdbname))
indexs = [val[1:] for val in upldf[(upldf['cya'] == 0)].index.tolist() if val[0] != 'P']
for x in range(0,len(indexs),50):
  i = x
  plmout = 'color c0, noes and resi '
  cmxout = 'color #{:}:'.format(mn)
  for j in range(50):
    plmout = plmout + str(indexs[i]) + '+'
    cmxout = cmxout + str(indexs[i]) + ','
    i=i+1
    if i== len(indexs): break
  outpml.write(plmout[:-1] + '\n')
  outcmx.write(cmxout[:-1] + ' c0 target ac\n')
for n in range(2,10,2):
  indexs = [val[1:] for val in upldf[(upldf['cya'] == n-1)].index.tolist()]
  indexs.extend([val[1:] for val in upldf[(upldf['cya'] == n)].index.tolist()])
  for x in range(0,len(indexs),50):
    i = x
    plmout = 'color c{:}, noes and resi '.format(str(n))
    cmxout = 'color #{:}:'.format(mn)
    for j in range(50):
      plmout = plmout + str(indexs[i]) + '+'
      cmxout = cmxout + str(indexs[i]) + ','
      i=i+1
      if i== len(indexs): break
    outpml.write(plmout[:-1] + '\n')
    outcmx.write(cmxout[:-1] + ' c{:} target ac\n'.format(str(n)))
indexs = [val[1:] for val in upldf[(upldf['cya'] >= 9)].index.tolist()]
for x in range(0,len(indexs),50):
  i = x
  plmout = 'color c10, noes and resi '
  cmxout = 'color #{:}:'.format(mn)
  for j in range(50):
    plmout = plmout + str(indexs[i]) + '+'
    cmxout = cmxout + str(indexs[i]) + ','
    i=i+1
    if i== len(indexs): break
  outpml.write(plmout[:-1] + '\n')
  outcmx.write(cmxout[:-1] + ' c10 target ac\n')
outcmx.write(sidechains[:-1].replace('#1',"#{:}".format(mn)) + '\n')
outcmx.write("show #{:}:thr,met,ala,leu,val,ile,phe,tyr\nhide #{:}@H*\ncolor byhetero\n".format(mn,mn))
outcmx.write('key c0:0 c2:2 c4:4 c6:6 c8:8 c10:10 fontsize 14 colorTreatment distinct numericLabelSpacing equal\nkey size 0.25000,0.03000\n')
outpml.write('show sticks, noes and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\nhide sticks, elem H\ncolor blue, elem N\ncolor gold, elem S\ncolor red, elem O\ncolor orange, elem P\ncolor white, elem H\nshow sticks, name N+H\n')
outpml.write(sidechains[:-1].replace(",","+").replace('#1:',"sticks, noes and resi ") + '\n')
outpml.close()
outcmx.close()

print('finished')

