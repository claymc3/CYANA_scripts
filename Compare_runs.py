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
'ALA':{'HA':'CA','CA':'HA','H':'N','N':'H','CB':'QB','QB':'CB'},
'CYS':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','H':'N','N':'H'},
'CYSS':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','H':'N','N':'H'},
'ASP':{'HA':'CA','HB2':'CB','CB':'HB2','HB3':'CB','CB':'HB3','QB':'CB','H':'N','N':'H'},
'GLU':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB3,HB2','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','H':'N','N':'H'},
'PHE':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB3,HB2','QB':'CB,HB2,HB3','QD':'CD1,CD2,HD1,HD2','QE':'CE1,CE2,HE1,HE2','HD1':'CD1','HE1':'CE1','HZ':'CZ','HE2':'CE2','HD2':'CD2','CD1':'HD1','CE1':'HE1','CZ': 'HZ','CE2':'HE2','CD2':'HD2','H':'N','N':'H'},
'GLY':{'QA':'CA,HA2,HA3','CA':'HA2,HA3','HA2':'CA','HA3':'CA','H':'N','N':'H'},
'HIS':{'HA':'CA','CA':'HA','CB':'HB2,HB3','HB2':'CB','HB3':'CB','QB':'CB,HB2,HB3','HD1':'ND1','HE2':'NE2','HD2':'CD2','HE1':'CE1','ND1':'HD1','NE2':'HE2','CD2':'HD2','CE1':'HE1','H':'N','N':'H'},
'HIST':{'HA':'CA','CA':'HA','CB':'HB2,HB3','HB2':'CB','HB3':'CB','QB':'CB,HB2,HB3','HD1':'ND1','HE2':'NE2','HD2':'CD2','HE1':'CE1','ND1':'HD1','NE2':'HE2','CD2':'HD2','CE1':'HE1','H':'N','N':'H'},
'HIS+':{'HA':'CA','CA':'HA','CB':'HB2,HB3','HB2':'CB','HB3':'CB','QB':'CB,HB2,HB3','HD1':'ND1','HE2':'NE2','HD2':'CD2','HE1':'CE1','ND1':'HD1','NE2':'HE2','CD2':'HD2','CE1':'HE1','H':'N','N':'H'},
'ILE':{'HA':'CA','HB':'CB','CA':'HA','CB':'HB','QG2':'CG2','QG1':'CG1,HG12,HG13','QD1':'CD1','CG2':'QG2','CG1':'QG1,HG12,HG13','CD1':'QD1','HG12':'CG1','HG13':'CG1','H':'N','N':'H'},
'LYS':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','HD2':'CD','HD3':'CD','CD':'HD2,HD3','QD':'CD,HD2,HD3','HE2':'CE','HE3':'CE','CE':'HE2,HE3','QE':'CE,HE2,HE3','QZ':'NZ','NZ':'QZ','H':'N','N':'H'},
'LEU':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','QD1':'CD1','CD1':'QD1','QD2':'CD2','CD2':'QD2','H':'N','N':'H'},
'MET':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','QE':'CE','CE':'QE','H':'N','N':'H'},
'ASN':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','QD2':'ND2,HD21,HD22','HD21':'ND2,QD2','HD22':'ND2,QD2','ND2':'HD21,HD22','H':'N','N':'H'},
'PRO':{'N':'N','HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','HD2':'CD','HD3':'CD','CD':'HD2,HD3','QD':'CD,HD2,HD3'},
'CPRO':{'N':'N','HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','HD2':'CD','HD3':'CD','CD':'HD2,HD3','QD':'CD,HD2,HD3'},
'GLN':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','HE21':'NE2','HE22':'NE2','QE2':'NE2,HE21,HE22','NE2':'HE21,HE22,QE2','H':'N','N':'H'},
'ARG':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG2':'CG','HG3':'CG','CG':'HG2,HG3','QG':'CG,HG2,HG3','HD2':'CD','HD3':'CD','CD':'HD2,HD3','QD':'CD,HD2,HD3','HE':'NE','NE':'HE','NH1':'HH11,HH12','HH11':'NH1','HH12':'NH1','QH1':'NH1,HH11,HH12','NH2':'HH21,HH22','HH21':'NH2','HH22':'NH2','QH2':'NH2,HH21,HH22','H':'N','N':'H'},
'SER':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG':'OG','OG':'HG','H':'N','N':'H'},
'SEP':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HG':'OG','OG':'HG','H':'N','N':'H'},
'THR':{'HA':'CA','CA':'HA','HB':'CB','CB':'HB','HG1':'OG1','OG1':'HG1','QG2':'CG2','CG2':'QG2','H':'N','N':'H'},
'TPO':{'HA':'CA','CA':'HA','HB':'CB','CB':'HB','HG1':'OG1','OG1':'HG1','QG2':'CG2','CG2':'QG2','H':'N','N':'H'},
'VAL':{'HA':'CA','CA':'HA','HB':'CB','CB':'HB','QG1':'CG1','CG1':'QG1','QG2':'CG2','CG2':'QG2','H':'N','N':'H'},
'TRP':{'HA':'CA','CA':'HA','HB2':'CB','HB3':'CB','CB':'HB2,HB3','QB':'CB,HB2,HB3','HD1':'CD1','CD1':'HD1','HE3':'CE3','HE1':'NE1','HZ3':'CZ3','HZ2':'CZ2','HH2':'CH2','CE3':'HE3','NE1':'HE1','CZ3':'HZ3','CZ2':'HZ2','CH2':'HH2','H':'N','N':'H'},
'TYR':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB3,HB2','QB':'CB,HB2,HB3','QD':'CD1,CD2,HD1,HD2','QE':'CE1,CE2,HE1,HE2','HD1':'CD1','HE1':'CE1','HZ':'CZ','HE2':'CE2','HD2':'CD2','CD1':'HD1','CE1':'HE1','CZ': 'HZ','CE2':'HE2','CD2':'HD2','HH':'OH','OH':'HH','H':'N','N':'H'},
'PTR':{'HA':'CA','HB2':'CB','HB3':'CB','CB':'HB3,HB2','QB':'CB,HB2,HB3','QD':'CD1,CD2,HD1,HD2','QE':'CE1,CE2,HE1,HE2','HD1':'CD1','HE1':'CE1','HZ':'CZ','HE2':'CE2','HD2':'CD2','CD1':'HD1','CE1':'HE1','CZ': 'HZ','CE2':'HE2','CD2':'HD2','HH':'OH','OH':'HH','H':'N','N':'H'}}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',
ConTypeDict = {'ALAH':'N', 'ALAHA':'Ali', 'ALAHB':'Methyl', 'ALAHB1':'Methyl', 'ALAHB2':'Methyl', 'ALAHB3':'Methyl', 'ALAQB':'Methyl', 'ALAC':'Other', 'ALACA':'Ali', 'ALACB':'Methyl', 'ALAN':'N', 'CYSSH':'N', 'CYSSHA':'Ali', 'CYSSHB2':'Ali', 'CYSSHB3':'Ali', 'CYSSQB':'Ali', 'CYSSC':'Other', 'CYSSCA':'Ali', 'CYSSCB':'Ali', 'CYSSN':'N', 'CYSH':'N', 'CYSHA':'Ali', 'CYSHB2':'Ali', 'CYSHB3':'Ali', 'CYSQB':'Ali', 'CYSHG':'Ali', 'CYSC':'Other', 'CYSCA':'Ali', 'CYSCB':'Ali', 'CYSN':'N', 'ASPH':'N', 'ASPHA':'Ali', 'ASPHB2':'Ali', 'ASPHB3':'Ali', 'ASPQB':'Ali', 'ASPHD2':'Other', 'ASPC':'Other', 'ASPCA':'Ali', 'ASPCB':'Ali', 'ASPCG':'Ali', 'ASPN':'N', 'GLUH':'N', 'GLUHA':'Ali', 'GLUHB2':'Ali', 'GLUHB3':'Ali', 'GLUQB':'Ali', 'GLUHE2':'Other', 'GLUHG2':'Ali', 'GLUHG3':'Ali', 'GLUQG':'Ali', 'GLUC':'Other', 'GLUCA':'Ali', 'GLUCB':'Ali', 'GLUCG':'Ali', 'GLUCD':'Other', 'GLUN':'N', 'PHEH':'N', 'PHEHA':'Aro', 'PHEHB2':'Aro', 'PHEHB3':'Aro', 'PHEQB':'Aro', 'PHEHD1':'Aro', 'PHEHD2':'Aro', 'PHEQD':'Aro', 'PHEHE1':'Aro', 'PHEHE2':'Aro', 'PHEQE':'Aro', 'PHEHZ':'Aro', 'PHEC':'Other', 'PHECA':'Aro', 'PHECB':'Aro', 'PHECD1':'Aro', 'PHECD2':'Aro', 'PHECE1':'Aro', 'PHECE2':'Aro', 'PHECG':'Aro', 'PHECZ':'Aro', 'PHEN':'N',  'GLYH':'N', 'GLYHA2':'Ali', 'GLYHA3':'Ali', 'GLYC':'Other', 'GLYCA':'Ali', 'GLYN':'N',  'HISH':'N', 'HISHA':'Ali', 'HISHB2':'Ali', 'HISHB3':'Ali', 'HISQB':'Ali', 'HISHD1':'N', 'HISHD2':'Aro', 'HISHE1':'Aro', 'HISC':'Other', 'HISCA':'Ali', 'HISCB':'Ali', 'HISCD2':'Aro', 'HISCE1':'Aro', 'HISCG':'Aro', 'HISN':'N', 'HISND1':'N', 'HISNE2':'N', 'HISTH':'N', 'HISTHA':'Ali', 'HISTHB2':'Ali', 'HISTHB3':'Ali', 'HISTQB':'Ali', 'HISTHD2':'Aro', 'HISTHE1':'Aro', 'HISTHE2':'N', 'HISTC':'Other', 'HISTCA':'Ali', 'HISTCB':'Ali', 'HISTCD2':'Aro', 'HISTCE1':'Aro', 'HISTCG':'Aro', 'HISTN':'N', 'HISTND1':'N', 'HISTNE2':'N', 'HIS+H':'N', 'HIS+HA':'Ali', 'HIS+HB2':'Ali', 'HIS+HB3':'Ali', 'HIS+QB':'Ali', 'HIS+HD1':'N', 'HIS+HD2':'Aro', 'HIS+HE1':'Aro', 'HIS+HE2':'N', 'HIS+C':'Other', 'HIS+CA':'Ali', 'HIS+CB':'Ali', 'HIS+CD2':'Aro', 'HIS+CE1':'Aro', 'HIS+CG':'Aro', 'HIS+N':'N', 'HIS+ND1':'N', 'HIS+NE2':'N', 'ILEH':'N', 'ILEHA':'Ali', 'ILEHB':'Ali', 'ILEHG12':'Ali', 'ILEHG13':'Ali', 'ILEQG1':'Ali', 'ILEHD1':'Methyl', 'ILEHD11':'Methyl', 'ILEHD12':'Methyl', 'ILEHD13':'Methyl', 'ILEQD1':'Methyl', 'ILEHG2':'Methyl', 'ILEHG21':'Methyl', 'ILEHG22':'Methyl', 'ILEHG23':'Methyl', 'ILEQG2':'Methyl', 'ILEC':'Other', 'ILECA':'Ali', 'ILECB':'Ali', 'ILECD1':'Methyl', 'ILECG1':'Ali', 'ILECG2':'Methyl', 'ILEN':'N', 'LYSH':'N', 'LYSHA':'Ali', 'LYSHB2':'Ali', 'LYSHB3':'Ali', 'LYSQB':'Ali', 'LYSHD2':'Ali', 'LYSHD3':'Ali', 'LYSQD':'Ali', 'LYSHE2':'Ali', 'LYSHE3':'Ali', 'LYSQE':'Ali', 'LYSHG2':'Ali', 'LYSHG3':'Ali', 'LYSQG':'Ali', 'LYSC':'Other', 'LYSCA':'Ali', 'LYSCB':'Ali', 'LYSCD':'Ali', 'LYSCE':'Ali', 'LYSCG':'Ali', 'LYSN':'N', 'LYSNZ':'N', 'LYSQZ':'N', 'LYSHZ':'N', 'LEUH':'N', 'LEUHA':'Ali', 'LEUHB2':'Ali', 'LEUHB3':'Ali', 'LEUQB':'Ali', 'LEUHG':'Ali', 'LEUHD1':'Methyl', 'LEUHD11':'Methyl', 'LEUHD12':'Methyl', 'LEUHD13':'Methyl', 'LEUQD1':'Methyl', 'LEUHD2':'Methyl', 'LEUHD21':'Methyl', 'LEUHD22':'Methyl', 'LEUHD23':'Methyl', 'LEUQD2':'Methyl', 'LEUC':'Other', 'LEUCA':'Ali', 'LEUCB':'Ali', 'LEUCG':'Ali', 'LEUCD1':'Methyl', 'LEUCD2':'Methyl', 'LEUN':'N', 'METH':'N', 'METHA':'Ali', 'METHB2':'Ali', 'METHB3':'Ali', 'METQB':'Ali', 'METHG2':'Ali', 'METHG3':'Ali', 'METQG':'Ali', 'METHE':'Methyl', 'METHE1':'Methyl', 'METHE2':'Methyl', 'METHE3':'Methyl', 'METQE':'Methyl', 'METC':'Other', 'METCA':'Ali', 'METCB':'Ali', 'METCE':'Methyl', 'METCG':'Ali', 'METN':'N', 'ASNH':'N', 'ASNHA':'Ali', 'ASNHB2':'Ali', 'ASNHB3':'Ali', 'ASNQB':'Ali', 'ASNHD21':'N', 'ASNHD22':'N', 'ASNQD':'N', 'ASNC':'Other', 'ASNCA':'Ali', 'ASNCB':'Ali', 'ASNCG':'Other', 'ASNN':'N', 'ASNND2':'N', 'PROHA':'Ali', 'PROHB2':'Ali', 'PROHB3':'Ali', 'PROQB':'Ali', 'PROHD2':'Ali', 'PROHD3':'Ali', 'PROQD':'Ali', 'PROHG2':'Ali', 'PROHG3':'Ali', 'PROQG':'Ali', 'PROC':'Other', 'PROCA':'Ali', 'PROCB':'Ali', 'PROCD':'Ali', 'PROCG':'Ali', 'PRON':'N', 'CPROHA':'Ali', 'CPROHB2':'Ali', 'CPROHB3':'Ali', 'CPROQB':'Ali', 'CPROHD2':'Ali', 'CPROHD3':'Ali', 'CPROQD':'Ali', 'CPROHG2':'Ali', 'CPROHG3':'Ali', 'CPROQG':'Ali', 'CPROC':'Other', 'CPROCA':'Ali', 'CPROCB':'Ali', 'CPROCD':'Ali', 'CPROCG':'Ali', 'CPRON':'N', 'GLNH':'N', 'GLNHA':'Ali', 'GLNHB2':'Ali', 'GLNHB3':'Ali', 'GLNQB':'Ali', 'GLNHE21':'N', 'GLNHE22':'N', 'GLNQE2':'N', 'GLNHG2':'Ali', 'GLNHG3':'Ali', 'GLNQG':'Ali', 'GLNC':'Other', 'GLNCA':'Ali', 'GLNCB':'Ali', 'GLNCD':'Other', 'GLNCG':'Ali', 'GLNN':'N', 'GLNNE2':'N', 'ARGH':'N', 'ARGHA':'Ali', 'ARGHB2':'Ali', 'ARGHB3':'Ali', 'ARGQB':'Ali', 'ARGHD2':'Ali', 'ARGHD3':'Ali', 'ARGQD':'Ali', 'ARGHG2':'Ali', 'ARGHG3':'Ali', 'ARGQG':'Ali', 'ARGHH11':'N', 'ARGHH12':'N', 'ARGQH1':'N', 'ARGHH21':'N', 'ARGHH22':'N', 'ARGQH2':'N', 'ARGC':'Other', 'ARGCA':'Ali', 'ARGCB':'Ali', 'ARGCD':'Ali', 'ARGCG':'Ali', 'ARGCZ':'Ali', 'ARGN':'N', 'ARGNE':'N', 'ARGNH1':'N', 'ARGNH2':'N', 'ARGHE':'N', 'SERH':'N', 'SERHA':'Ali', 'SERHB2':'Ali', 'SERHB3':'Ali', 'SERQB':'Ali', 'SERHG':'Other', 'SERC':'Other', 'SERCA':'Ali', 'SERCB':'Ali', 'SERN':'N', 'SEPH':'N', 'SEPHA':'Ali', 'SEPHB2':'Ali', 'SEPHB3':'Ali', 'SEPQB':'Ali', 'SEPHG':'Other', 'SEPC':'Other', 'SEPCA':'Ali', 'SEPCB':'Ali', 'SEPN':'N', 'THRH':'N', 'THRHA':'Ali', 'THRHB':'Ali', 'THRHG1':'Other', 'THRHG2':'Methyl', 'THRQG2':'Methyl', 'THRC':'Other', 'THRCA':'Ali', 'THRCB':'Ali', 'THRCG2':'Methyl', 'THRN':'N', 'TPOH':'N', 'TPOHA':'Ali', 'TPOHB':'Ali', 'TPOHG1':'Other', 'TPOHG2':'Methyl', 'TPOQG2':'Methyl', 'TPOC':'Other', 'TPOCA':'Ali', 'TPOCB':'Ali', 'TPOCG2':'Methyl', 'TPON':'N', 'VALH':'N', 'VALHA':'Ali', 'VALHB':'Ali', 'VALHG1':'Methyl', 'VALQG1':'Methyl', 'VALHG2':'Methyl', 'VALQG2':'Methyl', 'VALC':'Other', 'VALCA':'Ali', 'VALCB':'Ali', 'VALCG1':'Methyl', 'VALCG2':'Methyl', 'VALN':'N', 'TRPH':'N', 'TRPHA':'Ali', 'TRPHB2':'Ali', 'TRPHB3':'Ali', 'TRPQB':'Ali', 'TRPHD1':'Aro', 'TRPHE1':'N', 'TRPHE3':'Aro', 'TRPHH2':'Aro', 'TRPHZ2':'Aro', 'TRPHZ3':'Aro', 'TRPC':'Other', 'TRPCA':'Ali', 'TRPCB':'Ali', 'TRPCD1':'Aro', 'TRPCD2':'Aro', 'TRPCE2':'Aro', 'TRPCE3':'Aro', 'TRPCG':'Aro', 'TRPCH2':'Aro', 'TRPCZ2':'Aro', 'TRPCZ3':'Aro', 'TRPN':'N', 'TRPNE1':'N', 'TYRH':'N', 'TYRHA':'Aro', 'TYRHB2':'Aro', 'TYRHB3':'Aro', 'TYRQB':'Aro', 'TYRHD1':'Aro', 'TYRHD2':'Aro', 'TYRQD':'Aro', 'TYRHE1':'Aro', 'TYRHE2':'Aro', 'TYRQE':'Aro', 'TYRHH':'Other', 'TYRC':'Other', 'TYRCA':'Aro', 'TYRCB':'Aro', 'TYRCD1':'Aro', 'TYRCD2':'Aro', 'TYRCE1':'Aro', 'TYRCE2':'Aro', 'TYRCG':'Aro', 'TYRCZ':'Aro', 'TYRN':'N', 'PTRH':'N', 'PTRHA':'Aro', 'PTRHB2':'Aro', 'PTRHB3':'Aro', 'PTRQB':'Aro', 'PTRHD1':'Aro', 'PTRHD2':'Aro', 'PTRQD':'Aro', 'PTRHE1':'Aro', 'PTRHE2':'Aro', 'PTRQE':'Aro', 'PTRHH':'Other', 'PTRC':'Other', 'PTRCA':'Aro', 'PTRCB':'Aro', 'PTRCD1':'Aro', 'PTRCD2':'Aro', 'PTRCE1':'Aro', 'PTRCE2':'Aro', 'PTRCG':'Aro', 'PTRCZ':'Aro', 'PTRN':'N', 'N-N':'N_N', 'Aro-Aro':'Aro_Aro', 'Methyl-Methyl':'Methyl_Methyl', 'Ali-Ali':'Ali_Ali', 'Aro-Methyl':'Methyl_Aro', 'Methyl-Aro':'Methyl_Aro', 'Methyl-N':'N_Methyl', 'N-Methyl':'N_Methyl', 'Aro-N':'N_Aro', 'N-Aro':'N_Aro', 'Methyl-Ali':'Ali_Methyl', 'Ali-Methyl':'Ali_Methyl', 'Aro-Ali':'Ali_Aro', 'Ali-Aro':'Ali_Aro', 'Ali-N':'N_Ali', 'N-Ali':'N_Ali', 'Other-Other':'Other','Other-N':'Other','N-Other':'Other','Other-Ali':'Other','Ali-Other':'Other','Other-Methyl':'Other','Methyl-Other':'Other', 'Other-Aro':'Other', 'Aro-Other':'Other'}
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}
Ambiguous = {'ARGQB':'HB2,HB3','ARGQG':'HG2,HG3','ARGQG':'HG2,HG3','ARGQH1':'HH11,HH12','ARGQH2':'HH21,HH22','ASNQB':'HB2,HB3','ASNQD2':'HD21,HD22','ASPQB':'HB2,HB3','CYSQB':'HB2,HB3','CYSSQB':'HB2,HB3','GLNQB':'HB2,HB3','GLNQG':'HG2,HG3','GLNQE2':'HE21,HE22','GLUQB':'HB2,HB3','GLUQG':'HG2,HG3','GLYQA':'HA2,HA3','HISQB':'HB2,HB3','HISTQB':'HB2,HB3','ILEQG1':'HG12,HG13','LEUQB':'HB2,HB3','LEUQQD':'QD1,QD2','LYSQB':'HB2,HB3','LYSQG':'HG2,HG3','LYSQD':'HD2,HD3','LYSQE':'HE2,HE3','METQB':'HB2,HB3','METQG':'HG2,HG3','PHEQB':'HB2,HB3','PHEQG':'HG2,HG3','PHEQD':'HD2,HD3','PHEQE':'HE1,HE2','PROQB':'HB2,HB3','PROQG':'HG2,HG3','PROQD':'HD2,HD3','SERQB':'HB2,HB3','TRPQB':'HB2,HB3','TYRQB':'HB2,HB3','TYRQG':'HG2,HG3','TYRQD':'HD2,HD3','TYRQE':'HE1,HE2','PTRQB':'HB2,HB3','PTRQG':'HG2,HG3','PTRQD':'HD2,HD3','PTRQE':'HE1,HE2','VALQQG':'QG1,QG2'}

cwd = os.getcwd() + '/'

run1dir = sys.argv[1]
run2dir = sys.argv[2]

if run1dir[-1] != '/': run1dir = run1dir + '/'
if run2dir[-1] != '/': run2dir = run2dir + '/'
changelog = open(run2dir + '{:}_vs_{:}_Changes_log.txt'.format(sys.argv[2],sys.argv[1]),'w')
changelog.write('Reporting changes in {:} relative to {:}\n\n'.format(run2dir,run1dir))
fupl1 = run1dir + 'final.upl'
fupl2 = run2dir + 'final.upl'
hbupl1 = run1dir + 'hbond.upl'
hbupl2 = run2dir + 'hbond.upl'
hblol1 = run1dir + 'hbond.lol'
hblol2 = run2dir + 'hbond.lol'
dihe1 = run1dir + 'dihed.aco'
dihe2 = run2dir + 'dihed.aco'

# cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line and not re.match('^\s*#', line)][0].split()[2].split(',')
print('Reorting changes to input files of {:} relative to {:}'.format(run2dir,run1dir))

hbonds1,hbonds2,lol1,lol2 ={},{},{},{}
for ln in ['1','2']:
	hbdict = eval('hbonds{:}'.format(ln))
	loldict = eval('lol{:}'.format(ln))
	uplf = eval('hbupl{:}'.format(ln))
	lolf = eval('hblol{:}'.format(ln))
	for line in open(uplf).readlines():
		if line.strip() and not re.match(r'^\s*#', line):
			cns = line.split()
			upl = '{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])
			if upl in hbdict.keys():
				print('Duplicate hbond entry {:} in {:}'.format(upl,uplf))
			if upl not in hbdict.keys():
				hbdict[upl] = [cns[6]]
	for line in open(lolf).readlines():
		if line.strip() and not re.match(r'^\s*#', line):
			cns = line.split()
			lol = '{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])
			if lol in hbdict.keys():
				hbdict['{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])].append(cns[6])
			if lol not in hbdict.keys():
				print('missing hbond lol entry {:}'.format(upl))
			if lol in loldict.keys():
				print('Duplicate hbond entry {:} in {:}'.format(lol,lolf))
			if lol not in loldict.keys():
				loldict[lol] =[cns[6]]


rmhb, addhb = '',''
for hb in hbonds1.keys():
	if hb not in hbonds2.keys():
		rmhb = rmhb + '   {:}\n'.format(hb)
for hb in hbonds2.keys():
	if hb not in hbonds1.keys():
		addhb = addhb + '   {:}\n'.format(hb)
if len(rmhb) > 0: 
	rmhb = 'Removed Hbonds:\n' + rmhb + '\n'
	print(rmhb)
	changelog.write(rmhb)
if len(addhb) > 0: 
	addhb = 'Added Hbonds:\n' + addhb + '\n'
	changelog.write(addhb)
	print(addhb)
dihedrals1,dihedrals2 = {},{}
for ln in ['1','2']:
	angdict = eval('dihedrals{:}'.format(ln))
	dihedf = eval('dihe{:}'.format(ln))
	for line in open(dihedf).readlines():
		if line.strip() and not re.match(r'^\s*#', line):
			cns = line.split()
			# ang = '{:>5}  {:<4} {:<5}  {:>8} {:>8}'.format(cns[0],cns[1],cns[2],cns[3],cns[4])
			ang = '{:>5}  {:<4} {:<5}'.format(cns[0],cns[1],cns[2])
			if ang in angdict.keys():
				print('Duplicate dihedral entry {:} in {:}'.format(ang,dihedf))
			if ang not in angdict.keys():
				angdict[ang] = '{:6.1f}{:8.1f}'.format(float(cns[3]),float(cns[4]))
rmdih,addih,changedih = '','',''
for ang in dihedrals1.keys():
	if ang not in dihedrals2.keys():
		rmdih = rmdih + '   {:}  {:}\n'.format(ang,dihedrals1[ang])
for ang in dihedrals2.keys():
	if ang not in dihedrals1.keys():
		addih = addih  + '  {:}  {:}\n'.format(ang,dihedrals2[ang])
for ang in dihedrals1.keys():
	if ang in dihedrals2.keys():
		low1, up1 = dihedrals1[ang].split()
		low2, up2 = dihedrals2[ang].split()
		if low1 != low2 or up1 != up2:
			changedih = changedih + '   {:}  {:}  to  {:}\n'.format(ang,dihedrals1[ang],dihedrals2[ang])
if len(rmdih) > 0:
	rmdih = 'Removed Dihedrals:\n' + rmdih + '\n'
	changelog.write(rmdih)
	print(rmdih)
if len(addih) > 0:
	addih = 'Added Dihedrals:\n' + addih + '\n'
	changelog.write(addih)
	print(addih)
if len(changedih) > 0:
	changedih = 'Changed Dihedral Bounds:\n' + changedih + '\n'
	changelog.write(changedih)
	print(changedih)

calc = run2dir + '/CALC.cya'
manualcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
uplsf = [con for con in manualcons if 'upl' in con and 'hbond' not in con]
lolsf = [con for con in manualcons if 'lol' in con and 'hbond' not in con]
cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line and not re.match(r'^\s*#', line)][0].split()[2].split(',')

## Check to see if both directories have the same files
distcons = []
for uplf in uplsf:
	if os.path.exists(run1dir+ uplf):
		distcons.append(uplf)
	if not os.path.exists(run2dir+ uplf):
		changelog.write('{:} not used in {:}'.format(uplf,run1dir))
for lolf in lolsf:
	if os.path.exists(run1dir+ lolf):
		distcons.append(lolf)

for distcon in distcons:
	for ln in ['1','2']:
		exec("{:}{:} = {{}}".format(distcon.split('.')[0],ln))
		condict = eval("{:}{:}".format(distcon.split('.')[0],ln))
		confile = eval('run{:}dir'.format(ln)) + distcon 
		for line in open(confile).readlines():
			if line.strip() and not re.match(r'^\s*#', line):
				cns = line.split()
				upl = '{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])
				if upl in condict.keys():
					print('Duplicate entry {:} in {:}'.format(upl,confile))
				if upl not in condict.keys():
					condict['{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])] = float(cns[6])

for distcon in distcons:
	print(distcon)
	lostupl,newupl,changupl,chang2upl = '','','',''
	dict1 = eval(distcon.split('.')[0]+'1')
	dict2 = eval(distcon.split('.')[0]+'2')
	for upl in dict1.keys():
		if upl not in dict2.keys():
			try:
				atoms1 = replacements[upl.split()[1]][upl.split()[2]]
			except KeyError:
				atoms1 = upl.split()[2]
			try:
				atoms2 = replacements[upl.split()[4]][upl.split()[5]]
			except KeyError:
				atoms2 = upl.split()[5]
			badupl = ''
			for atom1 in atoms1.split(','):
				for atom2 in atoms2.split(','):
					rupl = '{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(upl.split()[0],upl.split()[1],atom1,upl.split()[4],upl.split()[5],atom2)
					if rupl not in dict2.keys():
						badupl = '   {:}  {:}\n'.format(upl,dict1[upl])
					if rupl in dict2.keys():
						changupl = changupl + '   {:}  {:} to {:} {:}\n'.format(upl,dict1[upl],rupl,dict2[rupl])
			if len(badupl) > 0: lostupl = lostupl +badupl
	for upl in dict2.keys():
		if upl not in dict1.keys():
			atoms1 = replacements[upl.split()[1]][upl.split()[2]]
			atoms2 = replacements[upl.split()[4]][upl.split()[5]]
			badupl = ''
			for atom1 in atoms1.split(','):
				for atom2 in atoms2.split(','):
					rupl = '{:>5}  {:<4}  {:<5}  {:>5}  {:<4}  {:<4}'.format(upl.split()[0],upl.split()[1],atom1,upl.split()[4],upl.split()[5],atom2)
					if rupl not in dict1.keys():
						badupl = '   {:}  {:}\n'.format(upl,dict2[upl])
					if rupl in dict1.keys():
						changupl = changupl + '   {:}  {:} to {:} {:}\n'.format(upl,dict2[upl],rupl,dict1[rupl])
				if len(badupl) > 0: newupl = newupl + badupl
	for upl in dict1.keys():
		if upl in dict2.keys():
			if float(dict1[upl]) != float(dict2[upl]):
				chang2upl = chang2upl + '   {:}  {:} to {:}\n'.format(upl,dict1[upl],dict2[upl])

	if len(lostupl) > 0:
		lostupl = 'Removed UPLs from {:}:\n'.format(distcon.split('.')[0]) + lostupl + '\n'
		print(lostupl)
		changelog.write(lostupl)
	if len(newupl) > 0:
		newupl = 'Added UPLs to {:}:\n'.format(distcon.split('.')[0]) + newupl + '\n'
		print(newupl)
		changelog.write(newupl)
	if len(changupl) > 0:
		changupl = 'Changed UPL Atoms in {:}:\n'.format(distcon.split('.')[0]) + changupl + '\n'
		print(changupl)
		changelog.write(changupl)
	if len(chang2upl) > 0:
		chang2upl = 'Changed UPL Value in {:}:\n'.format(distcon.split('.')[0]) + chang2upl + '\n'
		print(chang2upl)
		changelog.write(chang2upl)

#cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line and not re.match('^\s*#', line)][0].split()[2].split(',')


## Check if there are final.upl files available, if so compare them

if os.path.exists(fupl1) and os.path.exists(fupl2):
	finalupl1,finalupl2 ={},{}
	for ln in ['1','2']:
		upldict = eval('finalupl{:}'.format(ln))
		uplf = eval('fupl{:}'.format(ln))
		for line in open(uplf).readlines():
			if line.strip() and not re.match(r'^\s*#', line):
				cns = line.split()
				upl = '{:>} {:<4} {:<5} {:>} {:<4} {:<5}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])
				if upl not in upldict.keys():
					upldict['{:>} {:<4} {:<5} {:>} {:<4} {:<5}'.format(cns[0],cns[1],cns[2],cns[3],cns[4],cns[5])] = float(cns[6])
	lostupl,newupl,changupl = '','',''
	for upl in finalupl1.keys():
		if upl not in finalupl2.keys():
			lostupl = lostupl + '   {:}   {:6.2f}\n'.format(upl,finalupl1[upl])
	for upl in finalupl2.keys():
		if upl not in finalupl1.keys():
			newupl = newupl + '   {:}   {:6.2f}\n'.format(upl,finalupl2[upl])
	for upl in finalupl1.keys():
		if upl in finalupl2.keys():
			if abs(finalupl1[upl] - finalupl2[upl]) > 0.4:
				changupl = changupl + '   {:}  {:6.2f}  to  {:6.2f}\n'.format(upl,finalupl1[upl],finalupl2[upl])
	if len(lostupl) > 0:
		lostupl = '{:} {:} final UPLs not found in {:}:\n'.format(lostupl.count('\n'),run1dir,run2dir) + lostupl + '\n'
		print('{:} {:} final UPLs not found in {:}'.format(lostupl.count('\n'),run1dir,run2dir))
		changelog.write(lostupl)
	if len(newupl) > 0:
		newupl = '{:} {:} final UPLs not found in {:}:\n'.format(newupl.count('\n'),run2dir,run1dir) + newupl + '\n'
		print('{:} {:} final UPLs not found in {:}'.format(newupl.count('\n'),run2dir,run1dir))
		changelog.write(newupl)
	if len(changupl) > 0:
		changupl = '{:} UPL with difference > 0.4:\n'.format(changupl.count('\n')) + changupl + '\n'
		print('{:} UPL with difference > 0.4'.format(changupl.count('\n')))
		changelog.write(changupl)










