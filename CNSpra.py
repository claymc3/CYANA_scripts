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
import GetDistances as distance
from statistics import mean
# datetime object containing current date and time
now = datetime.now()
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y-%m-%d %H:%M")

Vnum = '1.1'
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB','CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB','CYSHB':'CB','CYSSHA':'CA','CYSSHB2':'CB','CYSSHB3':'CB','CYSSQB':'CB','CYSSHB':'CB','ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB','ASPHB':'CB','GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG','GLUHG':'CG','PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEHB':'CB','PHEQD':'CD1,CD2','PHEHD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2','GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA','GLYHA':'CA','HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1','HISTHA':'CA','HISTHB2':'CB','HISTHB3':'CB','HISTQB':'CB','HISTHB':'CB','HISTHD1':'ND1','HISTHE2':'NE2','HISTHD2':'CD2','HISTHE1':'CE1','HIS+HA':'CA','HIS+HB2':'CB','HIS+HB3':'CB','HIS+QB':'CB','HIS+HB':'CB','HIS+HD1':'ND1','HIS+HE2':'NE2','HIS+HD2':'CD2','HIS+HE1':'CE1','ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEHG1':'CG1','ILEQD1':'CD1','ILEHD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1','LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ','LYSHZ':'NZ','LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD1':'CD1','LEUHD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUHD2':'CD2','LEUQQD':'CD2,CD1','METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHB':'CB','METHG2':'CG','METHG3':'CG','METHG':'CG','METQG':'CG','METHG':'CG','METQE':'CE','METHE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE','ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD21':'ND2','ASNQD22':'ND2','ASNQD2':'ND2','ASNHD2':'ND2','PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD','PROHD':'CD','CPROHA':'CA','CPROHB2':'CB','CPROHB3':'CB','CPROQB':'CB','CPROHB':'CB','CPROHG2':'CG','CPROHG3':'CG','CPROQG':'CG','CPROHG':'CG','CPROHD2':'CD','CPROHD3':'CD','CPROQD':'CD','CPROHD':'CD','GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNHQ21':'NE2','GLNQE22':'NE2','GLNQE2':'NE2','GLNHE2':'NE2','ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2','ARGHH2':'NH2','SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG','SEPHA':'CA','SEPHB2':'CB','SEPHB3':'CB','SEPHG':'OG','THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2','THRHG2':'CG2','TPOHA':'CA','TPOHB':'CB','TPOHG1':'OG1','TPOHG21':'CG2','TPOHG22':'CG2','TPOHG23':'CG2','TPOQG2':'CG2','TPOHG2':'CG2','VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALHG2':'CG2','VALQQG':'CG1,CG2','TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2','TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRHB':'CB','TYRQD':'CD1,CD2','TYRHD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH','PTRHA':'CA','PTRHB2':'CB','PTRHB3':'CB','PTRQB':'CB','PTRHB':'CB','PTRQD':'CD1,CD2','PTRHD':'CD1,CD2','PTRQE':'CE1,CE2','PTRHE':'CE1,CE2','PTRHD1':'CD1','PTRHE1':'CE1','PTRHE2':'CE2','PTRHD2':'CD2','PTRHH':'OH','GLUHG1':'CG','ILEHG11':'CG','LYSHG1':'CG','LYSHD1':'CD','LYSHE1':'CE','METHG1':'CG','PROHG1':'CG','PROHD1':'CD','GLNHG1':'CG','ARGHG1':'CG','ARGHD1':'CD'}
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',
Pseudo2Prot = {'AQB':['HB1', 'HB2', 'HB3'],'RQB':['HB2', 'HB3'],'RQG':['HG2', 'HG3'],'RQD':['HD2', 'HD3'],'RQH1':['HH11', 'HH12'],'RQH2':['HH21', 'HH22'],'NQB':['HB2', 'HB3'],'NQD2':['HD21', 'HD22'],'DQB':['HB2', 'HB3'],'CQB':['HB2', 'HB3'],'QQB':['HB2', 'HB3'],'QQG':['HG2', 'HG3'],'QQE2':['HE21', 'HE22'],'EQB':['HB2', 'HB3'],'EQG':['HG2', 'HG3'],'GQA':['HA2', 'HA3'],'HQB':['HB2', 'HB3'],'IQG1':['HG11', 'HG12', 'HG13'],'IQG2':['HG21', 'HG22', 'HG23'],'IQD1':['HD11', 'HD12', 'HD13'],'LQB':['HB2', 'HB3'],'LQD1':['HD11', 'HD12', 'HD13'],'LQD2':['HD21', 'HD22', 'HD23'],'KQG':['HG2', 'HG3'],'KQD':['HD2', 'HD3'],'KQE':['HE2', 'HE3'],'KQZ':['HZ2', 'HZ3'],'MQB':['HB2', 'HB3'],'MQG':['HG2', 'HG3'],'MQE':['HE1', 'HE2', 'HE3'],'FQB':['HB2', 'HB3'],'FQD':['HD1', 'HD2'],'FQE':['HE1', 'HE2'],'PQB':['HB2', 'HB3'],'PQG':['HG2', 'HG3'],'PQD':['HD2', 'HD3'],'SQB':['HB2', 'HB3'],'TQG2':['HG21', 'HG22', 'HG23'],'WQB':['HB2', 'HB3'],'YQB':['HB2', 'HB3'],'YQD':['HD1', 'HD2'],'YQE':['HE1', 'HE2'],'VQB':['HB2', 'HB3'],'VQG1':['HG11', 'HG12', 'HG13'],'VQG2':['HG21', 'HG22', 'HG23'],'AHB':['HB1', 'HB2', 'HB3'],'RHB':['HB2', 'HB3'],'RHG':['HG2', 'HG3'],'RHD':['HD2', 'HD3'],'RHH1':['HH11', 'HH12'],'RHH2':['HH21', 'HH22'],'NHB':['HB2', 'HB3'],'NHD2':['HD21', 'HD22'],'DHB':['HB2', 'HB3'],'HHB':['HB2', 'HB3'],'QHB':['HB2', 'HB3'],'QHG':['HG2', 'HG3'],'QHE2':['HE21', 'HE22'],'EHB':['HB2', 'HB3'],'EHG':['HG2', 'HG3'],'GHA':['HA2', 'HA3'],'HHB':['HB2', 'HB3'],'IHG1':['HG12', 'HG13'],'IHG2':['HG21', 'HG22', 'HG23'],'IHD1':['HD11', 'HD12', 'HD13'],'LHB':['HB2', 'HB3'],'LHD1':['HD11', 'HD12', 'HD13'],'LHD2':['HD21', 'HD22', 'HD23'],'KHG':['HG2', 'HG3'],'KHD':['HD2', 'HD3'],'KHE':['HE2', 'HE3'],'KHZ':['HZ2', 'HZ3'],'MHB':['HB2', 'HB3'],'MHG':['HG2', 'HG3'],'MHE':['HE1', 'HE2', 'HE3'],'FHB':['HB2', 'HB3'],'FHD':['HD1', 'HD2'],'FHE':['HE1', 'HE2'],'PHB':['HB2', 'HB3'],'PHG':['HG2', 'HG3'],'PHD':['HD2', 'HD3'],'SHB':['HB2', 'HB3'],'THG2':['HG21', 'HG22', 'HG23'],'WHB':['HB2', 'HB3'],'YHB':['HB2', 'HB3'],'YHD':['HD1','HD2'],'YHE':['HE1', 'HE2'],'VHB':['HB2', 'HB3'],'VHG1':['HG11', 'HG12', 'HG13'],'VHG2':['HG21', 'HG22', 'HG23']}
ConTypeDict = {'ALAH':'N', 'ALAHA':'Ali', 'ALAHB':'Methyl', 'ALAHB1':'Methyl', 'ALAHB2':'Methyl', 'ALAHB3':'Methyl', 'ALAQB':'Methyl', 'ALAC':'Other', 'ALACA':'Ali', 'ALACB':'Methyl', 'ALAN':'N', 'CYSSH':'N', 'CYSSHA':'Ali', 'CYSSHB2':'Ali', 'CYSSHB3':'Ali', 'CYSSQB':'Ali', 'CYSSC':'Other', 'CYSSCA':'Ali', 'CYSSCB':'Ali', 'CYSSN':'N', 'CYSH':'N', 'CYSHA':'Ali', 'CYSHB2':'Ali', 'CYSHB3':'Ali', 'CYSQB':'Ali', 'CYSHG':'Ali', 'CYSC':'Other', 'CYSCA':'Ali', 'CYSCB':'Ali', 'CYSN':'N', 'ASPH':'N', 'ASPHA':'Ali', 'ASPHB2':'Ali', 'ASPHB3':'Ali', 'ASPQB':'Ali', 'ASPHD2':'Other', 'ASPC':'Other', 'ASPCA':'Ali', 'ASPCB':'Ali', 'ASPCG':'Ali', 'ASPN':'N', 'GLUH':'N', 'GLUHA':'Ali', 'GLUHB2':'Ali', 'GLUHB3':'Ali', 'GLUQB':'Ali', 'GLUHE2':'Other', 'GLUHG2':'Ali', 'GLUHG3':'Ali', 'GLUQG':'Ali', 'GLUC':'Other', 'GLUCA':'Ali', 'GLUCB':'Ali', 'GLUCG':'Ali', 'GLUCD':'Other', 'GLUN':'N', 'PHEH':'N', 'PHEHA':'Aro', 'PHEHB2':'Aro', 'PHEHB3':'Aro', 'PHEQB':'Aro', 'PHEHD1':'Aro', 'PHEHD2':'Aro', 'PHEQD':'Aro', 'PHEHE1':'Aro', 'PHEHE2':'Aro', 'PHEQE':'Aro', 'PHEHZ':'Aro', 'PHEC':'Other', 'PHECA':'Aro', 'PHECB':'Aro', 'PHECD1':'Aro', 'PHECD2':'Aro', 'PHECE1':'Aro', 'PHECE2':'Aro', 'PHECG':'Aro', 'PHECZ':'Aro', 'PHEN':'N',  'GLYH':'N', 'GLYHA2':'Ali', 'GLYHA3':'Ali', 'GLYC':'Other', 'GLYCA':'Ali', 'GLYN':'N',  'HISH':'N', 'HISHA':'Ali', 'HISHB2':'Ali', 'HISHB3':'Ali', 'HISQB':'Ali', 'HISHD1':'N', 'HISHD2':'Aro', 'HISHE1':'Aro', 'HISC':'Other', 'HISCA':'Ali', 'HISCB':'Ali', 'HISCD2':'Aro', 'HISCE1':'Aro', 'HISCG':'Aro', 'HISN':'N', 'HISND1':'N', 'HISNE2':'N', 'HISTH':'N', 'HISTHA':'Ali', 'HISTHB2':'Ali', 'HISTHB3':'Ali', 'HISTQB':'Ali', 'HISTHD2':'Aro', 'HISTHE1':'Aro', 'HISTHE2':'N', 'HISTC':'Other', 'HISTCA':'Ali', 'HISTCB':'Ali', 'HISTCD2':'Aro', 'HISTCE1':'Aro', 'HISTCG':'Aro', 'HISTN':'N', 'HISTND1':'N', 'HISTNE2':'N', 'HIS+H':'N', 'HIS+HA':'Ali', 'HIS+HB2':'Ali', 'HIS+HB3':'Ali', 'HIS+QB':'Ali', 'HIS+HD1':'N', 'HIS+HD2':'Aro', 'HIS+HE1':'Aro', 'HIS+HE2':'N', 'HIS+C':'Other', 'HIS+CA':'Ali', 'HIS+CB':'Ali', 'HIS+CD2':'Aro', 'HIS+CE1':'Aro', 'HIS+CG':'Aro', 'HIS+N':'N', 'HIS+ND1':'N', 'HIS+NE2':'N', 'ILEH':'N', 'ILEHA':'Ali', 'ILEHB':'Ali', 'ILEHG12':'Ali', 'ILEHG13':'Ali', 'ILEQG1':'Ali', 'ILEHD1':'Methyl', 'ILEHD11':'Methyl', 'ILEHD12':'Methyl', 'ILEHD13':'Methyl', 'ILEQD1':'Methyl', 'ILEHG2':'Methyl', 'ILEHG21':'Methyl', 'ILEHG22':'Methyl', 'ILEHG23':'Methyl', 'ILEQG2':'Methyl', 'ILEC':'Other', 'ILECA':'Ali', 'ILECB':'Ali', 'ILECD1':'Methyl', 'ILECG1':'Ali', 'ILECG2':'Methyl', 'ILEN':'N', 'LYSH':'N', 'LYSHA':'Ali', 'LYSHB2':'Ali', 'LYSHE1':'Ali', 'LYSHD1':'Ali', 'LYSHG1':'Ali', 'LYSHB1':'Ali', 'LYSHB3':'Ali', 'LYSQB':'Ali', 'LYSHD2':'Ali', 'LYSHD3':'Ali', 'LYSQD':'Ali', 'LYSHE2':'Ali', 'LYSHE3':'Ali', 'LYSQE':'Ali', 'LYSHG2':'Ali', 'LYSHG3':'Ali', 'LYSQG':'Ali', 'LYSC':'Other', 'LYSCA':'Ali', 'LYSCB':'Ali', 'LYSCD':'Ali', 'LYSCE':'Ali', 'LYSCG':'Ali', 'LYSN':'N', 'LYSNZ':'N', 'LYSQZ':'N', 'LYSHZ':'N', 'LEUH':'N', 'LEUHA':'Ali', 'LEUHB2':'Ali', 'LEUHB3':'Ali', 'LEUQB':'Ali', 'LEUHG':'Ali', 'LEUHD1':'Methyl', 'LEUHD11':'Methyl', 'LEUHD12':'Methyl', 'LEUHD13':'Methyl', 'LEUQD1':'Methyl', 'LEUHD2':'Methyl', 'LEUHD21':'Methyl', 'LEUHD22':'Methyl', 'LEUHD23':'Methyl', 'LEUQD2':'Methyl', 'LEUC':'Other', 'LEUCA':'Ali', 'LEUCB':'Ali', 'LEUCG':'Ali', 'LEUCD1':'Methyl', 'LEUCD2':'Methyl', 'LEUN':'N', 'METH':'N', 'METHA':'Ali', 'METHB2':'Ali', 'METHB3':'Ali', 'METQB':'Ali', 'METHG2':'Ali', 'METHG3':'Ali', 'METQG':'Ali', 'METHE':'Methyl', 'METHE1':'Methyl', 'METHE2':'Methyl', 'METHE3':'Methyl', 'METQE':'Methyl', 'METC':'Other', 'METCA':'Ali', 'METCB':'Ali', 'METCE':'Methyl', 'METCG':'Ali', 'METN':'N', 'ASNH':'N', 'ASNHA':'Ali', 'ASNHB2':'Ali', 'ASNHB3':'Ali', 'ASNQB':'Ali', 'ASNHD21':'N', 'ASNHD22':'N', 'ASNQD':'N', 'ASNC':'Other', 'ASNCA':'Ali', 'ASNCB':'Ali', 'ASNCG':'Other', 'ASNN':'N', 'ASNND2':'N', 'PROHA':'Ali', 'PROHB2':'Ali', 'PROHB3':'Ali', 'PROQB':'Ali', 'PROHD2':'Ali', 'PROHD3':'Ali', 'PROQD':'Ali', 'PROHG2':'Ali', 'PROHG3':'Ali', 'PROQG':'Ali', 'PROC':'Other', 'PROCA':'Ali', 'PROCB':'Ali', 'PROCD':'Ali', 'PROCG':'Ali', 'PRON':'N', 'CPROHA':'Ali', 'CPROHB2':'Ali', 'CPROHB3':'Ali', 'CPROQB':'Ali', 'CPROHD2':'Ali', 'CPROHD3':'Ali', 'CPROQD':'Ali', 'CPROHG2':'Ali', 'CPROHG3':'Ali', 'CPROQG':'Ali', 'CPROC':'Other', 'CPROCA':'Ali', 'CPROCB':'Ali', 'CPROCD':'Ali', 'CPROCG':'Ali', 'CPRON':'N', 'GLNH':'N', 'GLNHA':'Ali', 'GLNHB2':'Ali', 'GLNHB3':'Ali', 'GLNQB':'Ali', 'GLNHE21':'N', 'GLNHE22':'N', 'GLNQE2':'N', 'GLNHG2':'Ali', 'GLNHG3':'Ali', 'GLNQG':'Ali', 'GLNC':'Other', 'GLNCA':'Ali', 'GLNCB':'Ali', 'GLNCD':'Other', 'GLNCG':'Ali', 'GLNN':'N', 'GLNNE2':'N', 'ARGH':'N', 'ARGHA':'Ali', 'ARGHG1':'Ali', 'ARGHD1':'Ali', 'ARGHB1':'Ali', 'ARGHB2':'Ali', 'ARGHB3':'Ali', 'ARGQB':'Ali', 'ARGHD2':'Ali', 'ARGHD3':'Ali', 'ARGQD':'Ali', 'ARGHG2':'Ali', 'ARGHG3':'Ali', 'ARGQG':'Ali', 'ARGHH11':'N', 'ARGHH12':'N', 'ARGQH1':'N', 'ARGHH21':'N', 'ARGHH22':'N', 'ARGQH2':'N', 'ARGC':'Other', 'ARGCA':'Ali', 'ARGCB':'Ali', 'ARGCD':'Ali', 'ARGCG':'Ali', 'ARGCZ':'Ali', 'ARGN':'N', 'ARGNE':'N', 'ARGNH1':'N', 'ARGNH2':'N', 'ARGHE':'N', 'SERH':'N', 'SERHA':'Ali', 'SERHB2':'Ali', 'SERHB3':'Ali', 'SERQB':'Ali', 'SERHG':'Other', 'SERC':'Other', 'SERCA':'Ali', 'SERCB':'Ali', 'SERN':'N', 'SEPH':'N', 'SEPHA':'Ali', 'SEPHB2':'Ali', 'SEPHB3':'Ali', 'SEPQB':'Ali', 'SEPHG':'Other', 'SEPC':'Other', 'SEPCA':'Ali', 'SEPCB':'Ali', 'SEPN':'N', 'THRH':'N', 'THRHA':'Ali', 'THRHB':'Ali', 'THRHG1':'Other', 'THRHG2':'Methyl', 'THRQG2':'Methyl', 'THRC':'Other', 'THRCA':'Ali', 'THRCB':'Ali', 'THRCG2':'Methyl', 'THRN':'N', 'TPOH':'N', 'TPOHA':'Ali', 'TPOHB':'Ali', 'TPOHG1':'Other', 'TPOHG2':'Methyl', 'TPOQG2':'Methyl', 'TPOC':'Other', 'TPOCA':'Ali', 'TPOCB':'Ali', 'TPOCG2':'Methyl', 'TPON':'N', 'VALH':'N', 'VALHA':'Ali', 'VALHB':'Ali', 'VALHG1':'Methyl', 'VALQG1':'Methyl', 'VALHG2':'Methyl', 'VALQG2':'Methyl', 'VALC':'Other', 'VALCA':'Ali', 'VALCB':'Ali', 'VALCG1':'Methyl', 'VALCG2':'Methyl', 'VALN':'N', 'TRPH':'N', 'TRPHA':'Ali', 'TRPHB2':'Ali', 'TRPHB3':'Ali', 'TRPQB':'Ali', 'TRPHD1':'Aro', 'TRPHE1':'N', 'TRPHE3':'Aro', 'TRPHH2':'Aro', 'TRPHZ2':'Aro', 'TRPHZ3':'Aro', 'TRPC':'Other', 'TRPCA':'Ali', 'TRPCB':'Ali', 'TRPCD1':'Aro', 'TRPCD2':'Aro', 'TRPCE2':'Aro', 'TRPCE3':'Aro', 'TRPCG':'Aro', 'TRPCH2':'Aro', 'TRPCZ2':'Aro', 'TRPCZ3':'Aro', 'TRPN':'N', 'TRPNE1':'N', 'TYRH':'N', 'TYRHA':'Aro', 'TYRHB2':'Aro', 'TYRHB3':'Aro', 'TYRQB':'Aro', 'TYRHD1':'Aro', 'TYRHD2':'Aro', 'TYRQD':'Aro', 'TYRHE1':'Aro', 'TYRHE2':'Aro', 'TYRQE':'Aro', 'TYRHH':'Other', 'TYRC':'Other', 'TYRCA':'Aro', 'TYRCB':'Aro', 'TYRCD1':'Aro', 'TYRCD2':'Aro', 'TYRCE1':'Aro', 'TYRCE2':'Aro', 'TYRCG':'Aro', 'TYRCZ':'Aro', 'TYRN':'N', 'PTRH':'N', 'PTRHA':'Aro', 'PTRHB2':'Aro', 'PTRHB3':'Aro', 'PTRQB':'Aro', 'PTRHD1':'Aro', 'PTRHD2':'Aro', 'PTRQD':'Aro', 'PTRHE1':'Aro', 'PTRHE2':'Aro', 'PTRQE':'Aro', 'PTRHH':'Other', 'PTRC':'Other', 'PTRCA':'Aro', 'PTRCB':'Aro', 'PTRCD1':'Aro', 'PTRCD2':'Aro', 'PTRCE1':'Aro', 'PTRCE2':'Aro', 'PTRCG':'Aro', 'PTRCZ':'Aro', 'PTRN':'N', 'N-N':'N_N', 'Aro-Aro':'Aro_Aro', 'Methyl-Methyl':'Methyl_Methyl', 'Ali-Ali':'Ali_Ali', 'Aro-Methyl':'Methyl_Aro', 'Methyl-Aro':'Methyl_Aro', 'Methyl-N':'N_Methyl', 'N-Methyl':'N_Methyl', 'Aro-N':'N_Aro', 'N-Aro':'N_Aro', 'Methyl-Ali':'Ali_Methyl', 'Ali-Methyl':'Ali_Methyl', 'Aro-Ali':'Ali_Aro', 'Ali-Aro':'Ali_Aro', 'Ali-N':'N_Ali', 'N-Ali':'N_Ali', 'Other-Other':'Other','Other-N':'Other','N-Other':'Other','Other-Ali':'Other','Ali-Other':'Other','Other-Methyl':'Other','Methyl-Other':'Other', 'Other-Aro':'Other', 'Aro-Other':'Other'}
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}
Ambiguous = {'ARGQB':'HB2,HB3','ARGQG':'HG2,HG3','ARGQG':'HG2,HG3','ARGQH1':'HH11,HH12','ARGQH2':'HH21,HH22','ASNQB':'HB2,HB3','ASNQD2':'HD21,HD22','ASPQB':'HB2,HB3','CYSQB':'HB2,HB3','CYSSQB':'HB2,HB3','GLNQB':'HB2,HB3','GLNQG':'HG2,HG3','GLNQE2':'HE21,HE22','GLUQB':'HB2,HB3','GLUQG':'HG2,HG3','GLYQA':'HA2,HA3','HISQB':'HB2,HB3','HISTQB':'HB2,HB3','ILEQG1':'HG12,HG13','LEUQB':'HB2,HB3','LEUQQD':'QD1,QD2','LYSQB':'HB2,HB3','LYSQG':'HG2,HG3','LYSQD':'HD2,HD3','LYSQE':'HE2,HE3','METQB':'HB2,HB3','METQG':'HG2,HG3','PHEQB':'HB2,HB3','PHEQG':'HG2,HG3','PHEQD':'HD2,HD3','PHEQE':'HE1,HE2','PROQB':'HB2,HB3','PROQG':'HG2,HG3','PROQD':'HD2,HD3','SERQB':'HB2,HB3','TRPQB':'HB2,HB3','TYRQB':'HB2,HB3','TYRQG':'HG2,HG3','TYRQD':'HD2,HD3','TYRQE':'HE1,HE2','PTRQB':'HB2,HB3','PTRQG':'HG2,HG3','PTRQD':'HD2,HD3','PTRQE':'HE1,HE2','VALQQG':'QG1,QG2'}

if len(sys.argv)==1:
  print('''

Usage: 
    cnspra [name]

Required Input:

    name        assuming your running script in CNS/refinePDB/ directory 
                name is the name you supplied to the cnsprep and cnspost
                scripts 

OutPut:
    post_cns_analysis/
        name.cxc
        name.pml
        name_summary.txt
        name_overview.pdf
        pdb_Phi.csv, pdb_Psi.csv, pdb_Chi1.csv, pdb_Chi2.csv 
        pseudobonds/ (Groups of Distance restraints)
            upl_constraint_type.pb
            upl_viol_upls.pb: UPLs violated in 10 or more structures 
            hbond.pb
''')
  exit()
colors = ['royalblue','forest','yellowgreen', 'darkorange','purple','lightseagreen ','darkkhaki','peru','saddlebrown','mediumpurple','blue','purple']
ConectionTypes = ['N_N','N_Methyl','N_Aro','Methyl_Methyl','Methyl_Aro','Aro_Aro','N_Ali','Ali_Ali','Ali_Aro','Ali_Methyl','Other','Intramolecular']
## Set up the directory locations to retreve information assuming that useres are running in the CNS/refinedPDB directory 
## cwd will be the CNS/refinedPDB; topdir willl be the calculation directory where the CNS directory is listed, 
name = sys.argv[1]
cwd = os.getcwd() + '/'
topdir = ''
for dirn in cwd.split('/'):
  if re.match('CNS',dirn): 
    cnsdir = dirn
    break
  else:
    topdir = topdir + dirn + '/'
outdir = topdir + 'post_cns_analysis/'
in_pdb = '{:}_cya.pdb'.format(name)
noetbl = open('../{:}_noe.tbl'.format(name)).readlines()
hbondtbl = open('../{:}_hbond.tbl'.format(name)).readlines()
dihetbl = open('../{:}_dihe.tbl'.format(name)).readlines()
pdbname = '{:}_cya'.format(name)
outname = '{:}_cns'.format(name)
init = '../../init.cya'
seq = [line.strip().split() for line in open(topdir + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0] and line.split()[0] in AAA_dict.keys()]
Seqdict, num2AAA, ASequence = {}, {}, []
for resn,resi in seq:
  if resn in AAA_dict.keys():
    Seqdict[resi] = AAA_dict[resn] + resi
    num2AAA[resi] = resn
    ASequence.append(AAA_dict[resn] + resi)
upldf = pd.DataFrame(index = ASequence, columns=['cns','viol < 0.3','viol > 0.3'])
upldf['cns'] = np.zeros(len(ASequence))
upldf['viol < 0.3'] = np.zeros(len(ASequence))
upldf['viol > 0.3'] = np.zeros(len(ASequence))

## Check for the output directory if it does not exist make it
if not os.path.exists(outdir):
  os.makedirs(outdir)
if not os.path.exists(outdir +'pseudobonds/'):
  os.makedirs(outdir +'pseudobonds/')
checkcons = open(outdir + outname + '_summary.txt','w')

shortsum = open(outdir + 'Short_stats.txt','w')
checkcons.write('## Generated using CNSpra_{:} on {:} \n'.format(Vnum,dt_string))
shortsum.write('## Generated using CNSpra_{:} on {:} \n'.format(Vnum,dt_string))

outpml = open('{:}{:}.pml'.format(outdir,outname),'w')
outpml.write('load ../{:}/refinedPDB/{:}_cya.pdb\n'.format(cnsdir,name))
outpml.write('set dash_gap, 0.05\n')
outpml.write('set_color royalblue = [65,105,225]\nset_color forest = [34,139,34]\nset_color yellowgreen = [154,205,50]\nset_color darkorange = [255,140,0]\nset_color purple = [128,0,128]\nset_color lightseagreen = [32,178,170]\nset_color darkkhaki = [189,183,107]\nset_color peru = [205,133,63]\nset_color saddlebrown = [139,69,19]\nset_color gold = [255,215,0]\nset_color navy = [0,0,128]\nset_color darkturquoise = [0,206,209]\nset_color pink = [255,192,203]\nset_color cyan = [0,255,255]\nset_color paleturquoise = [175,238,238]\nset_color lightsalmon = [255,160,122]\nset_color khaki = [240,230,140]\nset_color yellowgreen = [154,205,50]\nset_color thistle = [216,191,216]\nset_color aquamarine = [127,255,212]\nset_color plum = [221,160,221]\nset_color lightpink = [255,182,193]\nset_color mediumvioletred = [199,21,133]\nset_color firebrick = [178,34,34]\nset_color lightcoral = [240,128,128]\nset_color deeppink = [255,20,147]\nset_color hotpink = [255,105,180]\nset_color mediumpurple = [147,112,219]\nset_color navy = [0,0,128]\nset_color cornflowerblue = [100,149,237]\n')
outpml.write('color gray60, all\n')
outpml.write('show sticks, {:} and resn THR+MET+ALA+LEU+VAL+ILE+PHE+TYR\n hide sticks, elem H\nhide sticks, name N+C\n'.format(pdbname))
outpml.write('color paleturquoise, {:} and resn ILE\ncolor lightsalmon, {:} and resn LEU\ncolor khaki, {:} and resn VAL\ncolor yellowgreen, {:} and resn ALA\ncolor thistle, {:} and resn MET\ncolor aquamarine, {:} and resn THR\ncolor lightpink, {:} and resn TYR\ncolor plum, {:} and resn PHE\n'.format(pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname,pdbname))
outpml.write('color gold, elem S\ncolor red, elem O\ncolor blue, elem N\n')
outpml.write('split_states ' + pdbname + '\n')
for y in range(2,21,1):
  outpml.write('align {:}_{:04d}, {:}_0001\n'.format(pdbname,y, pdbname))
outcmx = open('{:}{:}.cxc'.format(outdir,name),'w')
outcmx.write('open ../{:}/refinedPDB/{:}_cya.pdb\n'.format(cnsdir,name))
outcmx.write('color #1 gray(150)\n')
outcmx.write('match #1.2-20 to #1.1\n')
outcmx.write('color #1:ile paleturquoise target a\ncolor #1:leu lightsalmon  target a\ncolor #1:val khaki target a\ncolor #1:ala yellowgreen  target a\ncolor #1:met thistle target a\ncolor #1:thr aquamarine target a\ncolor #1:phe plum target a\ncolor #1:tyr lightpink target a\n')
outcmx.write('show #1:thr,met,ala,leu,val,ile,phe,tyr\nname meyfside #1:thr,met,ala,leu,val,ile,phe,tyr\ncartoon suppress false\n')
outcmx.write('label #1.1 text "{0.label_one_letter_code}{0.number}{0.insertion_code}"\n''label ontop false\n')
outcmx.write('ui tool show "Side View"\n#ui mousemode right distance\n')
## ---------------------------------------------------------------------------
## Make pseudo bond and gorup statements for poor, long, and short upl 
## entries from designated upl file 

for conect in ConectionTypes:
    exec("{:}_pb = []".format(conect))
    exec("group{:} = 'group {:}, '".format(conect, conect))

viol1list, viol2list = [],[]
viol1pbout = open(outdir +'pseudobonds/dist_viol_<_0.3.pb','w')
viol1pbout.write("; halfbond = false\n; color = mediumvioletred\n; radius = 0.1\n; dashes = 0\n")
viol2pbout = open(outdir +'pseudobonds/dist_viol_>_0.3.pb','w')
viol2pbout.write("; halfbond = false\n; color = firebrick\n; radius = 0.1\n; dashes = 0\n")
viol1cons, viol2cons = 'group dis_viol_0.1-0.3, ', 'group dist_viol_0.3, '
violdict, dihedviol = {}, {}

## ---------------------------------------------------------------------------
## Make a dictionary for each model key = A#-atom value = [x,y,z] 

Starts,Ends = [], []
for mnum in range(1,21,1):
  start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
  exec('Coor' + str(mnum) + ' = {}')
  Starts.append(start)
  Ends.append(start-1)
Ends.append(len(open(in_pdb).readlines()))
Ends = Ends[1:]
pdb = open(in_pdb).readlines()
n = 0
for (start,end) in zip(Starts,Ends):
  n+=1
  Coor = eval('Coor' + str(n))
  for x in range(start,end,1):
    line = pdb[x]
    if line[0:4] == "ATOM" or line[0:4] == 'HETA':
      if line[17:20].strip() in AAA_dict.keys():
        index = '{:}{:}-{:}'.format(AAA_dict[line[17:20].strip()],line[22:26].strip(),line[12:16].strip())
        Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

## ---------------------------------------------------------------------------
##  Read in name_noe.tbl file and extract violations > 0.1 that occur in 5 or 
##  more structures 

DistViolDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'max viol','average','num viol'])
i, im, v = 1, 0, 0
sidelist = []
for noe in noetbl:
  if not re.match('^\s+\!',noe):
    noe = noe.replace('HN)','H').replace('HN ','H ').replace('HB1','HB3').replace('(','').replace(')','')
    if re.search("\\*",noe.split()[5]): 
      noea1 = noe.split()[5].replace('H','Q').replace('*','')
      noe = noe.replace(noe.split()[5],noea1)
    if re.search("\\*",noe.split()[10]): 
      noea2 = noe.split()[10].replace('H','Q').replace('*','')
      noe = noe.replace(noe.split()[10],noea2)
    cns = noe.split()
    if num2AAA[cns[2]]+cns[5] in ConTypeDict.keys():ct1 = ConTypeDict[num2AAA[cns[2]]+cns[5]]
    if num2AAA[cns[2]]+cns[5] not in ConTypeDict.keys(): ct1 = 'Other'
    if num2AAA[cns[7]]+cns[10] in ConTypeDict.keys():ct2 = ConTypeDict[num2AAA[cns[7]]+cns[10]]
    if num2AAA[cns[7]]+cns[10] not in ConTypeDict.keys(): ct2 = 'Other'
    ctype = ConTypeDict["{:}-{:}".format(ct1,ct2)]
    pblist = eval(ctype + '_pb')
    atom1 = cns[5]
    atom2 = cns[10]
    if num2AAA[cns[2]]+cns[5] in replacements.keys():
      atom1 = atom1.replace(cns[5], replacements[num2AAA[cns[2]]+cns[5]])
    if num2AAA[cns[7]]+cns[10] in replacements.keys():
      atom2 = atom2.replace(cns[10], replacements[num2AAA[cns[7]]+cns[10]])
    upldf.loc[Seqdict[cns[2]],'cns'] = upldf.loc[Seqdict[cns[2]],'cns'] + 1
    upldf.loc[Seqdict[cns[7]],'cns'] = upldf.loc[Seqdict[cns[7]],'cns'] + 1
    ## Traslate amides H to N to search input upls but don't mess up the connections drawn in pymol/chimera
    if atom1 == 'H': atom1 = 'N'
    if atom2 == 'H': atom2 = 'N'
    if (float(cns[2]) > 1000 and float(cns[7]) < 1000) or (float(cns[2]) < 1000 and float(cns[7]) > 1000):
      im+=1
      ctype = 'Intramolecular'
      pblist = eval('Intramolecular_pb')
      outpml.write('distance intramol{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(im), pdbname, cns[2], atom1, pdbname, cns[7], atom2))
      exec('groupIntramolecular = groupIntramolecular+ "intramo{:} "'.format(im))
      pblist.append('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[2], atom1, cns[7],atom2))
    ## Make sure all connections to side chains are shown
    if (num2AAA[cns[2]] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[5] not in ['N','H']) and cns[2] not in sidelist:
      sidelist.append(cns[2])
    if (num2AAA[cns[7]] not in ['ALA','LEU','VAL','MET','ILE','THR','TYR','PHE'] and cns[10] not in ['N','H']) and cns[7] not in sidelist:
      sidelist.append(cns[7])
    i+=1
    outpml.write('distance UPL{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(i), pdbname, cns[2], atom1, pdbname, cns[7], atom2))
    pblist.append('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[2], atom1, cns[7],atom2))
    exec('group' + ctype + '=' + 'group' + ctype + '+ "UPL{:} "'.format(i))
    dmin = np.round(float(cns[11]) - float(cns[12]),3)
    dmax = np.round(float(cns[11]) + float(cns[13]),3)
    cons1 = '{:}-{:}-{:}-{:}'.format(Seqdict[cns[2]],cns[5],Seqdict[cns[7]],cns[10])
    Coord = eval('Coor1')
    matom1 = distance.find_protons('{:}-{:}'.format(Seqdict[cns[2]],cns[5]), Coord)
    matom2 = distance.find_protons('{:}-{:}'.format(Seqdict[cns[7]],cns[10]), Coord)
    dist = []
    for mnum in range(1,21,1):
      Coor = eval('Coor' + str(mnum))
      d = distance.getDistance(matom1, matom2, Coor)
      if np.round(d - dmax,3) >= 0.100:
        DistViolDF.loc[cons1,mnum] = np.round(d - dmax,2)
        dist.append(np.round(d - dmax,2))
      if np.round(dmin - d,3) >= 0.100:
        DistViolDF.loc[cons1,mnum] = np.round(d - dmin,2)
    if len(dist) > 5:
      DistViolDF.loc[cons1,'max viol'] = max(dist)
      DistViolDF.loc[cons1,'average'] = np.round(mean(dist),2)
      DistViolDF.loc[cons1,'num viol'] = len(dist)
      v+=1
      outpml.write('distance viol{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(v), pdbname, cns[2], atom1, pdbname, cns[7], atom2))
      violdict[cons1] = ""
      if max(dist) < 0.30:
        upldf.loc[Seqdict[cns[2]],'viol < 0.3'] = upldf.loc[Seqdict[cns[2]],'viol < 0.3'] + 1
        upldf.loc[Seqdict[cns[7]],'viol < 0.3'] = upldf.loc[Seqdict[cns[7]],'viol < 0.3'] + 1
        viol1cons = viol1cons + '+ "viol{:} "'.format(v)
        viol1pbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[2], atom1, cns[7],atom2))
        viol1list.append("{:<20}  {:3.2f}A violated in {:2} by {:3.2f}\u00C5\n".format(cons1,dmax, len(dist),np.round(mean(dist),2)))
      if max(dist) >= 0.30:
        upldf.loc[Seqdict[cns[2]],'viol > 0.3'] = upldf.loc[Seqdict[cns[2]],'viol > 0.3'] + 1
        upldf.loc[Seqdict[cns[7]],'viol > 0.3'] = upldf.loc[Seqdict[cns[7]],'viol > 0.3'] + 1
        viol2cons = viol2cons + '+ "viol{:} "'.format(v)
        viol2pbout.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[2], atom1, cns[7],atom2))
        viol2list.append("{:<20}  {:3.2f}A violated in {:2} by {:3.2f}\u00C5\n".format(cons1,dmax, len(dist),np.round(mean(dist),2)))
DistViolDF.to_csv(outdir + name + '_dist_viols_CNS.csv')
if im != 0:
  print('found {:} intermolecuar connections'.format(im))
DVioltext = '{:3.0f} Violated Distance < 0.3\u00C5\n{:3.0f} Violated Distance > 0.3\u00C5\n'.format(len(
  viol1list),len(viol2list))
print(DVioltext)

mn = 1
for x in range(len(ConectionTypes)):
  pbs = eval('{:}_pb'.format(ConectionTypes[x]))
  if len(pbs) > 1:
    mn+=1
    pbout = open('{:}pseudobonds/{:}_{:}.pb'.format(outdir, outname, ConectionTypes[x]),'w')
    pbout.write("; halfbond = false\n; color = " + colors[x] + "\n; radius = 0.1\n; dashes = 0\n")
    pbout.writelines(pbs)
    outcmx.write('open pseudobonds/{:}_{:}.pb\n'.format(outname, ConectionTypes[x]))
    groupstr = eval('group' + ConectionTypes[x])
    outpml.write(groupstr + '\n')
    outpml.write('color {:}, {:}\n'.format(colors[x],ConectionTypes[x]))
for (group,gname, color) in [('viol1cons','dist_viol_<_0.3','mediumvioletred'),('viol2cons','dist_viol_>_0.3','firebrick')]:
  mn+=1
  outcmx.write('open pseudobonds/{:}.pb\n'.format(gname))
  grpstr = eval(group)
  outpml.write(grpstr + '\n')
  outpml.write('color {:}, {:}\n'.format(color, gname))
selhbond = 'name hbond #1.1:'
hbonsl = []
hbond = open(outdir +'pseudobonds/hbond.pb','w')
hbond.write("; halfbond = false\n; color = pink\n; radius = 0.2\n; dashes = 10\n")
hbgroupline = 'group hbond , '
h = 1
mn+=1
for line in hbondtbl:
  cns = line.replace('HN','H').replace('(','').replace(')','').split()
  if not re.match('^\s*\!',line):
    if line.strip() and "#" not in cns[2]:
      if (cns[2],cns[7],cns[10]) not in hbonsl:
        h+=1 
        hbonsl.append((cns[2],cns[7],cns[10].replace('N','H')))
        cons = '{:}-{:}-{:}-{:}'.format(Seqdict[cns[2]],cns[5],Seqdict[cns[7]],cns[10])
        # if cons in Upperdict.keys():
        #   hbond.write('#1.1:{:}@{:} #1.1:{:}@{:} hotpink\n'.format(cns[2], cns[5], cns[7],cns[10]))
        # else:
        hbond.write('#1.1:{:}@{:} #1.1:{:}@{:}\n'.format(cns[2], cns[5], cns[7],cns[10]))
        outpml.write('distance hbond{:}, {:}_0001 and resi {:} and name {:}, {:}_0001 and resi {:} and name {:}\n'.format(str(h), pdbname, cns[2], cns[5].replace('H','N'), pdbname, cns[7], cns[10].replace('H','N')))
        hbgroupline = hbgroupline + 'hbond' + str(h) + ' '
        if cns[2] not in selhbond:
          selhbond = selhbond +'{:},'.format(cns[2])
        if cns[7] not in selhbond:
          selhbond = selhbond +'{:},'.format(cns[7])
        if cns[5] not in ['O','H','N','HN']: sidelist.append(cns[2])
        if cns[10] not in ['O','H','N','HN']: sidelist.append(cns[7])
hbond.close()
sidechains = 'show #1:'
for res in sidelist:
  sidechains = sidechains + res + ','
outcmx.write(sidechains[:-1] + '\n')
outpml.write('show sticks, resi '+ sidechains.split(':')[-1][:-1].replace(',','+')+'\n')
outcmx.write('hide H\n''show #1.1@H,N target a\n')
outpml.write(hbgroupline + '\n')
outpml.write('color pink, hbond\n')
outpml.write('show sticks,name N+O+C+H and resi '+ selhbond.split(':')[-1][:-1].replace(',','+')+'\n')
selhbond = selhbond[:-1] + '@O,N\nshow hbond target a\n'
outcmx.write('open pseudobonds/' + 'hbond.pb\n')
outcmx.write(selhbond)
## ---------------------------------------------------------------------------
## Manage the dihedral angles side of things 
phicount, psicount,chi1count,chi2count,othercount,total,vphicount, vpsicount,vchi1count,vchi2count,vothercount = 0,0,0,0,0,0,0,0,0,0,0
phipsidict,chidict,plotdict = {}, {}, {}
phiaco, chiaco, phiviol, chiviol, diheviols = [], [], [], [], []
AngleViolDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'max viol','average','num viol'])
for x in range(len(dihetbl)):
  line = dihetbl[x]
  if re.match('^\s*assign', line):
    if line.split()[2] not in Seqdict.keys():
      print('Bad dihedral entry\n{:}{:}'.format(line,dihetbl[x+1]))
      pass
    else:
      line2 = dihetbl[x+1].replace('HG21)','CD1 )').replace('CG2','CG1')
      total+=1
      mang = float(line2.split()[13])
      flex = float(line2.split()[14])
      a1 = '{:}-{:}'.format(Seqdict[line.split()[2]],line.split()[5])
      a2 = '{:}-{:}'.format(Seqdict[line.split()[8]],line.split()[11])
      a3 = '{:}-{:}'.format(Seqdict[line2.split()[1]],line2.split()[4])
      a4 = '{:}-{:}'.format(Seqdict[line2.split()[7]],line2.split()[10])
      resn = line.split()[8]
      if line.split()[5] == 'C' and line.split()[11] == 'N' and line2.split()[4] == 'CA' and line2.split()[10] == 'C': angle = 'PHI'
      if line.split()[5] == 'N' and line.split()[11] == 'CA' and line2.split()[4] == 'C' and line2.split()[10] == 'N': angle = 'PSI'
      if line.split()[5] == 'N' and line.split()[11] == 'CA' and line2.split()[4] == 'CB': angle = 'CHI1'
      if line.split()[5] == 'CA' and line.split()[11] == 'CB': angle = 'CHI2'
      try:
        exec('{:}count = {:}count + 1'.format(angle.lower(),angle.lower()))
      except NameError:
        othercount+=1
      if 'CHI' in angle and mang < 0: mang = mang + 360
      ang_min = np.round(float(mang) - float(line2.split()[14]),3)
      ang_min_orig = np.round(float(mang) - float(line2.split()[14]),3)
      ang_max = np.round(float(mang) + float(line2.split()[14]),3)
      if 'P' in angle: 
        angdict = phipsidict
        if resn not in phiaco:
          phiaco.append(resn)
      if 'CHI' in angle:
        angdict = chidict
        if ang_min < 0:
          ang_min = ang_min + 360.0
          ang_max = ang_max + 360.0
        if resn not in chiaco:
          chiaco.append(resn)
      outline = r"$\{:}$  {:} - {:}".format(angle.lower(), ang_min, ang_max)
      plotdict[Seqdict[resn] + angle] = [ang_min, ang_max]
      if Seqdict[resn] + anlge  not in angdict.keys():
        angdict[Seqdict[resn] + angle] = [[outline,'black']]
      else: 
        angdict[Seqdict[resn] + angle].append([outline,'black'])
      if ang_min <= -180.0:
        ang_min = ang_min + 360.0
        if ang_max < 0:
          ang_max = ang_max + 360.0
      # if ang_max > 180:

      angles = []
      for mnum in range(1,21,1):
        Coords = eval('Coor' + str(mnum))
        ang = np.round(Dihed.calcDihedrals(Coords[a1],Coords[a2],Coords[a3],Coords[a4]),1)
        if 'P' in angle and ang_min_orig < -180.0: ang = ang + 360.0
        if 'P' in angle and ang_max > 180.0 and ang < 0.0: ang = ang + 360.0
        # if ang < 0: ang = ang + 360.0
        if 'CHI' in angle and ang < 0.0: ang = ang + 360.0
        # if 'P' in angle and ang_max > 180:  ang_min = ang_max - 360.0 
        if np.round(ang_min - ang,1) > 5.0 :
          print('{:} max {:} min {:}'.format(num2AAA[resn]+resn+angle,ang_max,ang_min))
          print('min {:} val {:} min {:} diff {:}'.format(num2AAA[resn]+angle, ang, ang_min,np.round(ang_min - ang,1)))
          angles.append(np.round(ang_min - ang,1))
          AngleViolDF.loc[Seqdict[resn]+angle,mnum] = np.round(ang_min - ang,1)
        if np.round(ang - ang_max,1) > 5.0 : 
          angles.append(np.round(ang - ang_max,1))
          print('{:} max {:} min {:}'.format(num2AAA[resn]+resn+angle,ang_max,ang_min))
          print('max {:} val {:} max {:} diff {:}'.format(num2AAA[resn]+angle, ang, ang_max,np.round(ang - ang_max,1)))
          AngleViolDF.loc[Seqdict[resn]+angle,mnum] = np.round(ang - ang_max,1)
      if len(angles) > 5:
        diheviols.append("{:<6} {:<4s} Violated in {:2d} by {:3.1f}\n".format(Seqdict[resn],angle,len(angles),np.round(mean(angles),1)))
        AngleViolDF.loc[Seqdict[resn]+angle,'max viol'] = max(angles)
        AngleViolDF.loc[Seqdict[resn]+angle,'average'] = np.round(mean(angles),1)
        AngleViolDF.loc[Seqdict[resn]+angle,'num viol'] = len(angles)
        dihedviol[num2AAA[resn]+angle] = r'$\{:}$ viol in {:} by {:}'.format(angle.lower(), len(angles), np.round(mean(angles),1))
        try:
          exec('v{:}count = v{:}count + 1'.format(angle.lower(),angle.lower()))
        except NameError:
          vothercount+=1
        if angle == 'PHI' or angle == 'PSI':
          if resn not in phiviol:
            phiviol.append(resn)
        if 'CHI' in angle:
          if resn not in chiviol:
            chiviol.append(resn)
AngleViolDF.to_csv(outdir + name + '_dihed_viols_CNS.csv')
angle_text = "Total of {:} dihedral restraints:\n       input viol\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n{:<6} {:^5} {:^4}\n\n".format(total, 'Phi', phicount, vphicount, 'Psi', psicount, vpsicount , 'Chi1', chi1count ,vchi1count, 'Chi2', chi2count, vchi2count)
print(angle_text[:-2])
checkcons.write(angle_text)
shortsum.write(angle_text)

## ---------------------------------------------------------------------------
## Run the GetDihed.py to determine phi, psi, chi1 and chi2 and plot them
## for all 20 structures
## ---------------------------------------------------------------------------
text = [['UPL','#9acd32'],['Violated UPL < 0.3\u00C5','#ffa500'],['Violated UPL > 0.3\u00C5','#db7093']]
print('Extracting dihedrals')
DAramalist, DArotalist = Dihed.extract(in_pdb, ASequence, outdir, upldf, phipsidict, chidict, plotdict,dihedviol,text)
print('finished plotting dihedrals')
armn = mn+1
routmn = mn+2
mn+=2
outcmx.write('open ../{:}/refinedPDB/{:} maxModels 1\nrename #{:} angle_restraints\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(cnsdir,in_pdb,armn,armn,armn))
outcmx.write('open ../{:}/refinedPDB/{:} maxModels 1\nrename #{:} rama_outliers\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(cnsdir,in_pdb,routmn,routmn,routmn))
ramalist, rotalist = [],[]
for line in DAramalist:ramalist.append(line.split()[0][1:])
for line in DArotalist:rotalist.append(line.split()[0][1:])
anglesout = [[phiaco, "phipsi","purple",armn], [chiaco,"chi1chi2","navy",armn], 
[phiviol,"viol_phipsi","mediumpurple",armn], [chiviol,"chi1chi2","cornflowerblue",armn],
[ramalist,"dissallowed_phipsi","mediumvioletred",routmn], [rotalist,"dissallowed_chi1chi2","mediumvioletred",routmn]]
for listn, aname, color, modle in anglesout:
  outpml.write('create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\n'.format(aname,pdbname,aname,aname))
  rline = ""
  if len(listn) > 1:
    for resn in listn:
      rline = rline + "{:},".format(resn)
    rline = rline[:-1]
    if 'phi' in aname:
      outcmx.write("name {:} #{:}:{:}\ncolor {:} {:} target c\n".format(aname,modle,rline,aname,color))
      outpml.write("create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\ncolor {:}, {:} and resn {:}".format(aname,pdbname,aname,aname,color,aname,rline.replace(',','+')))
    if 'chi' in aname:
      outcmx.write("name {:} #{:}:{:}\nshow {:} target a\ncolor {:} {:} target a\nhide #{:}@H*,N,O target a\ncolor byhetero target a\n".format(aname,modle,rline,aname,aname,color,modle))
      outpml.write("create {:}, {:}_0001\ncolor gray60,{:}\nhide sticks, {:}\nshow sticks, {:} and resn {:}\ncolor {:}, {:} and resn {:}".format(aname,pdbname,aname,aname,aname,rline.replace(',','+'),color,aname,rline.replace(',','+')))

## ---------------------------------------------------------------------------
## Write things out the the summary file
## ---------------------------------------------------------------------------
AVioltext = '{:3.0f} Violated Dihedral Restraints\n{:3.0f} Residues with Disallowed Phi/Psi \n{:3.0f} Residues with Disallowed Chi1/Chi2 \n\n'.format(len(diheviols),len(DAramalist),len(DArotalist))
print(AVioltext)
checkcons.write(DVioltext)
checkcons.write(AVioltext)

checkcons.write('### {:3.0f}  Violated Distance Restraints  < 0.3 ###\n'.format(len(viol1list)))
# violpeaks = sorted(violpeaks, key = lambda x: (x.split()[10],x.split()[8]))
viol1list = sorted(viol1list, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
for viol in viol1list:
  # if viol in assigndict.keys():
    # checkcons.write('{:}  {:3.2f}A ({:}): {:}'.format(viol,float(upldict[viol]),len(assigndict[viol]), violdict[viol]))
    # checkcons.writelines(assigndict[viol])
  checkcons.write(viol)
  checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Poor/Low Support constraints to the summary file
viol2list = sorted(viol2list, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))
checkcons.write('### {:3.0f}  Violated Distance Restraints  > 0.3 \u00C5 ###\n'.format(len(viol2list)))
for viol in viol2list:
  # if con in assigndict.keys():
  #   checkcons.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
  #   checkcons.writelines(assigndict[con])
  checkcons.write(viol)
  checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Violated Dihedral Restraints to the summary file
checkcons.write('### {:3.0f} Violated Dihedral Restraints ###\n'.format(len(diheviols)))
diheviols = sorted(diheviols, key = lambda x: (x.split()[0][1:], x.split()[1]))
for con in diheviols:
  checkcons.writelines(con)
  checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Residues with Disallowed Phi/Psi to the summary file
checkcons.write('### {:3.0f} Disallowed Phi/Psi Dihedral Restraints ###\n'.format(len(DAramalist)))
for con in DAramalist:
  checkcons.writelines(con)
  checkcons.write('\n')
checkcons.write('\n\n')
#### Write out Residues with Disallowed Chi1/Chi2 to the summary file
checkcons.write('### {:3.0f} Disallowed Chi1/Chi2 Dihedral ###\n'.format(len(DArotalist)))
for con in DArotalist:
  checkcons.writelines(con)
  checkcons.write('\n')
checkcons.write('\n\n')
checkcons.close()
# exit()
# for line in open(fovw).readlines():
#   if line.strip():
#     if line.split()[0] == 'Ave':
#       shortsum.write("Target Function {:}\n".format(line.split()[1]))
#     if line.split()[0] == 'Average':
#       shortsum.write(line.strip()[8:52] + '\n')
# shortsum.close()
## ---------------------------------------------------------------------------
### Creating Model coloring residues based on the number of NOE restraints 
mn+=1
outcmx.write('open ../{:}/refinedPDB/{:}_cya.pdb maxModels 1\nrename #{:} noes\nhide #{:} target a\ncolor #{:} gray(150)\n'.format(cnsdir,name,mn,mn,mn))
outcmx.write('color name c0 rgb(255,205,0)\ncolor name c2 rgb(156,217,59)\ncolor name c4 rgb(52,182,121)\ncolor name c6 rgb(42,117,142)\ncolor name c8 rgb(59,81,139)\ncolor name c10 rgb(20,64,110)\n')
outpml.write('set_color c0 = [255,205,0]\nset_color c2  = [156,217,59]\nset_color c4  = [52,182,121]\nset_color c6  = [42,117,142]\nset_color c8  = [59,81,139]\nset_color c10  = [20,64,110]\n')
outpml.write('create noes, {:}_0001\ncolor gray60,phi-psi\nhide sticks, noes\n'.format(pdbname))
indexs = [val[1:] for val in upldf[(upldf['cns'] == 0)].index.tolist() if val[0] != 'P']
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
  indexs = [val[1:] for val in upldf[(upldf['cns'] == n-1)].index.tolist()]
  indexs.extend([val[1:] for val in upldf[(upldf['cns'] == n)].index.tolist()])
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
indexs = [val[1:] for val in upldf[(upldf['cns'] >= 9)].index.tolist()]
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

