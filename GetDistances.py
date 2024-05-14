'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
#
# ------------------------------------------------------------------------------

'''
import os
import sys
import pandas as pd
import numpy as np
import re
import math as m
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

Pseudo2Prot = {'AQB':['HB1', 'HB2', 'HB3'],'RQB':['HB2', 'HB3'],'RQG':['HG2', 'HG3'],'RQD':['HD2', 'HD3'],'RQH1':['HH11', 'HH12'],'RQH2':['HH21', 'HH22'],'NQB':['HB2', 'HB3'],'NQD2':['HD21', 'HD22'],'DQB':['HB2', 'HB3'],'CQB':['HB2', 'HB3'],'QQB':['HB2', 'HB3'],'QQG':['HG2', 'HG3'],'QQE2':['HE21', 'HE22'],'EQB':['HB2', 'HB3'],'EQG':['HG2', 'HG3'],'GQA':['HA2', 'HA3'],'HQB':['HB2', 'HB3'],'IQG1':['HG11', 'HG12', 'HG13'],'IQG2':['HG21', 'HG22', 'HG23'],'IQD1':['HD11', 'HD12', 'HD13'],'LQB':['HB2', 'HB3'],'LQD1':['HD11', 'HD12', 'HD13'],'LQD2':['HD21', 'HD22', 'HD23'],'KQG':['HG2', 'HG3'],'KQD':['HD2', 'HD3'],'KQE':['HE2', 'HE3'],'KQZ':['HZ2', 'HZ3'],'MQB':['HB2', 'HB3'],'MQG':['HG2', 'HG3'],'MQE':['HE1', 'HE2', 'HE3'],'FQB':['HB2', 'HB3'],'FQD':['HD1', 'HD2'],'FQE':['HE1', 'HE2'],'PQB':['HB2', 'HB3'],'PQG':['HG2', 'HG3'],'PQD':['HD2', 'HD3'],'SQB':['HB2', 'HB3'],'TQG2':['HG21', 'HG22', 'HG23'],'WQB':['HB2', 'HB3'],'YQB':['HB2', 'HB3'],'YQD':['HD1', 'HD2'],'YQE':['HE1', 'HE2'],'VQB':['HB2', 'HB3'],'VQG1':['HG11', 'HG12', 'HG13'],'VQG2':['HG21', 'HG22', 'HG23'],'AHB':['HB1', 'HB2', 'HB3'],'RHB':['HB2', 'HB3'],'RHG':['HG2', 'HG3'],'RHD':['HD2', 'HD3'],'RHH1':['HH11', 'HH12'],'RHH2':['HH21', 'HH22'],'NHB':['HB2', 'HB3'],'NHD2':['HD21', 'HD22'],'DHB':['HB2', 'HB3'],'HHB':['HB2', 'HB3'],'QHB':['HB2', 'HB3'],'QHG':['HG2', 'HG3'],'QHE2':['HE21', 'HE22'],'EHB':['HB2', 'HB3'],'EHG':['HG2', 'HG3'],'GHA':['HA2', 'HA3'],'HHB':['HB2', 'HB3'],'IHG1':['HG12', 'HG13'],'IHG2':['HG21', 'HG22', 'HG23'],'IHD1':['HD11', 'HD12', 'HD13'],'LHB':['HB2', 'HB3'],'LHD1':['HD11', 'HD12', 'HD13'],'LHD2':['HD21', 'HD22', 'HD23'],'KHG':['HG2', 'HG3'],'KHD':['HD2', 'HD3'],'KHE':['HE2', 'HE3'],'KHZ':['HZ2', 'HZ3'],'MHB':['HB2', 'HB3'],'MHG':['HG2', 'HG3'],'MHE':['HE1', 'HE2', 'HE3'],'FHB':['HB2', 'HB3'],'FHD':['HD1', 'HD2'],'FHE':['HE1', 'HE2'],'PHB':['HB2', 'HB3'],'PHG':['HG2', 'HG3'],'PHD':['HD2', 'HD3'],'SHB':['HB2', 'HB3'],'THG2':['HG21', 'HG22', 'HG23'],'WHB':['HB2', 'HB3'],'YHB':['HB2', 'HB3'],'YHD':['HD1','HD2'],'YHE':['HE1', 'HE2'],'VHB':['HB2', 'HB3'],'VHG1':['HG11', 'HG12', 'HG13'],'VHG2':['HG21', 'HG22', 'HG23']}
Prot2Heavy = {
'AHA':'CA','AQB':'CB','AHB':'CB','AHB1':'CB','AHB2':'CB','AHB3':'CB','AH':'N',
'CHA':'CA','CHB2':'CB','CHB3':'CB','CQB':'CB','CH':'N',
'DHA':'CA','DHB2':'CB','DHB3':'CB','DQB':'CB','DH':'N',
'EHA':'CA','EHB2':'CB','EHB3':'CB','EQB':'CB','EHG2':'CG','EHG3':'CG','EQG':'CG','EH':'N',
'FHA':'CA','FHB2':'CB','FHB3':'CB','FQB':'CB','FQD':'CD1,CD2','FQE':'CE1,CE2','FHD1':'CD1','FHE1':'CE1','FHZ':'CZ','FHE2':'CE2','FHD2':'CD2','FH':'N',
'GHA2':'CA','GHA3':'CA','GQA':'CA','GH':'N',
'HHA':'CA','HHB2':'CB','HHB3':'CB','HQB':'CB','HHD1':'ND1','HHE2':'NE2','HHD2':'CD2','HHE1':'CE1','HH':'N',
'IHA':'CA','IHB':'CB','IQG2':'CG2','IHG2':'CG2','IHG21':'CG2','IHG22':'CG2','IHG23':'CG2','IHG12':'CG1','IHG13':'CG1','IHG1':'CG1','IQG1':'CG1','IHD1':'CD1','IQD1':'CD1','IHD11':'CD1','IHD12':'CD1','IHD13':'CD1','IH':'N',
'KHA':'CA','KHB2':'CB','KHB3':'CB','KQB':'CB','KHG2':'CG','KHG3':'CG','KHD2':'CD ','KHD3':'CD ','KQD':'CD ','KHE2':'CE','KHE3':'CE','KQE':'CE','KHE':'CE','KHZ1':'NZ','KHZ2':'NZ','KHZ3':'NZ','KQZ':'NZ','KH':'N',
'LHA':'CA','LHB2':'CB','LHB3':'CB','LQB':'CB','LHG':'CG','LHD11':'CD1','LHD12':'CD1','LHD13':'CD1','LQD1':'CD1','LHD1':'CD1','LHD21':'CD2','LHD22':'CD2','LHD23':'CD2','LQD2':'CD2','LHD2':'CD2','LH':'N',
'MHA':'CA','MHB2':'CB','MHB3':'CB','MQB':'CB','MHG2':'CG','MHG3':'CG','MQG':'CG','MHE':'CE','MQE':'CE','MHE1':'CE','MHE2':'CE','MHE3':'CE','MH':'N',
'NHA':'CA','NHB2':'CB','NHB3':'CB','NQB':'CB','NHD21':'ND2','NHD22':'ND2','NQD2':'ND2','NH':'N',
'PHA':'CA','PHB2':'CB','PHB3':'CB','PQB':'CB','PHG2':'CG','PHG3':'CG','PQG':'CG','PHD2':'CD','PHD3':'CD','PQD':'CD',
'QHA':'CA','QHB2':'CB','QHB3':'CB','QQB':'CB','QHG2':'CG','QHG3':'CG','QQG':'CG','QHE21':'NE2','QHE22':'NE2','QQE2':'NE2','QH':'N',
'RHA':'CA','RHB2':'CB','RHB3':'CB','RQB':'CB','RHG2':'CG','RHG3':'CG','RQG':'CG','RHD2':'CD','RHD3':'CD','RQD':'CD','RHE':'NE','RHH11':'NH1','RHH12':'NH1','RQH1':'NH1','RHH21':'NH2','RHH22':'NH2','RQH2':'NH2','RH':'N',
'SHA':'CA','SHB2':'CB','SHB3':'CB','SHG':'OG','SH':'N',
'THA':'CA','THB':'CB','THG1':'OG1','THG21':'CG2','THG22':'CG2','THG23':'CG2','TQG2':'CG2','THG2':'CG2','TH':'N',
'VHA':'CA','VHB':'CB','VHG11':'CG1','VHG12':'CG1','VHG13':'CG1','VQG1':'CG1','VHG1':'CG1','VHG21':'CG2','VHG22':'CG2','VHG23':'CG2','VQG2':'CG2','VHG2':'CG2','VH':'N',
'WHA':'CA','WHB2':'CB','WHB3':'CB','WQB':'CB','WHD1':'CD1','WHE3':'CE3','WHE1':'NE1','WHZ3':'CZ3','WHZ2':'CZ2','WHH2':'CH2','WH':'N',
'YHA':'CA','YHB2':'CB','YHB3':'CB','YQB':'CB','YQD':'CD1,CD2','YQE':'CE1,CE2','YHD1':'CD1','YHE1':'CE1','YHE2':'CE2','YHD2':'CD2','YHH':'OH','YH':'N',}

AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS','I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

# cwd = '/Users/mclay1/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_70/'
# # cwd = '/Users/mclay1/FGFR3_structure/cyana_23/'

# in_pdb = cwd + 'CNS/refinedPDB/FGFR2_r70_cya.pdb'
# pdb_name = in_pdb.split('/')[-1].replace('.pdb','')
# fupl = cwd + 'final.upl'
# calc = cwd + 'CALC.cya'
# init = cwd + 'init.cya'


# prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line and not re.match('^\s*#', line)][0].split()[2].split(',')

def find_protons(inatom, PDBdict):
  atype = '{:}{:}'.format(inatom.split('-')[0][0], inatom.split('-')[-1])
  resid, name = inatom.split('-')[0],inatom.split('-')[-1]
  atoms = []
  if atype in Pseudo2Prot.keys():
    for atomid in ['{:}-{:}'.format(resid,atom) for atom in Pseudo2Prot[atype]]:
      if atomid in PDBdict.keys():
        atoms.append(atomid)
  if atype not in Pseudo2Prot.keys():
    if inatom in PDBdict.keys():
      atoms.append(inatom)
  return atoms
# ------------------------------------------------------------------------------
# return atoms_list entry containing heavy atom for any pseudo atom entry 
# to allow examination of the effect of sum over r^6 averaging 
# ------------------------------------------------------------------------------
def find_heavy(inatom, PDBdict):
  atype = '{:}{:}'.format(inatom.split('-')[0][0], inatom.split('-')[-1])
  resid, name = inatom.split('-')[0],inatom.split('-')[-1]
  atoms = []
  if atype in Prot2Heavy.keys():
    atomid = '{:}-{:}'.format(resid,Prot2Heavy[atype])
    if atomid in PDBdict.keys():
      atoms.append(atomid)
  return atoms
# ------------------------------------------------------------------------------
def getDistance(donor, acceptor,PDBdict):
  d = 0.0
  dist = []
  for a1 in donor:
    for a2 in acceptor:
      (x1,y1,z1) = PDBdict[a1]
      (x2,y2,z2) = PDBdict[a2]
      d = d + np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**-6
      dist.append(d)
  if d != 0.0: 
    reff = np.round(d**(-1/6),2)
    #reff = round(m.pow(d/len(dist),-1/6),1)
  else: reff = 0.0
  return reff

def examin(in_pdb, ADpairs,Assignments,finalupl):
  pdb_name = in_pdb.split('/')[-1].replace('.pdb','')
  num2seq = {}
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
    # print('Reading coordinates for model {:d}'.format(n))
    Coor = eval('Coor' + str(n))
    for x in range(start,end,1):
      line = pdb[x]
      if line[0:4] == "ATOM" or line[0:4] == 'HETA':
        if line[17:20].strip() in AAA_dict.keys():
          index = '{:}{:}-{:}'.format(AAA_dict[line[17:20].strip()],line[22:26].strip(),line[12:16].strip())
          Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
  DistancesDF = pd.DataFrame(columns=Assignments,index=Assignments)
  FilteredDF = pd.DataFrame(columns=Assignments,index=Assignments)
  # from itertools import combinations
  # ADpairs = [ '{:}-{:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
  for connection in ADpairs:
    inatom1 = '{:}-{:}'.format(connection.split('-')[0],connection.split('-')[1])
    inatom2 = '{:}-{:}'.format(connection.split('-')[2],connection.split('-')[3])
    Coord = eval('Coor1')
    atom1 = find_protons(inatom1, Coord)
    heavy1 = find_heavy(inatom1, Coord)
    atom2 = find_protons(inatom2, Coord)
    heavy2 = find_heavy(inatom2, Coord)
    dist = []
    # print(heavy1,heavy2)
    if heavy1 != heavy2:
      for mnum in range(1,21,1):
        Coor = eval('Coor' + str(mnum))
        d = getDistance(atom1, atom2, Coor)
        # if getDistance(atom1, atom2, Coor) <= 7.5:
        dist.append(getDistance(heavy1, heavy2, Coor))
      if len(dist) >= 5:
        DistancesDF.loc[inatom1,inatom2] = "{:3.2f} +/- {:0.2f}".format(np.mean(dist),np.std(dist))
        DistancesDF.loc[inatom2,inatom1] = "{:3.2f} +/- {:0.2f}".format(np.mean(dist),np.std(dist))
        if np.mean(dist) <= 7.5:
          FilteredDF.loc[inatom1,inatom2] = "{:3.2f} +/- {:0.2f}".format(np.mean(dist),np.std(dist))
          FilteredDF.loc[inatom2,inatom1] = "{:3.2f} +/- {:0.2f}".format(np.mean(dist),np.std(dist))
          # if np.std(dist) >= 0.5:
          #   print("{:} {:} {:3.2f} +/- {:0.2f}".format(inatom1,inatom2,np.mean(dist),np.std(dist)))
  cyana,observed,cyana1 = [],[],[]
  for line in open(finalupl).readlines():
    if '#SUP' not in line: ## exclude ambiguous restraints
      pass 
    else:
      cns = line.split()
      atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
      atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
      upl = float(cns[6])
      cyana.append(upl)
      cyana1.append(upl-1.0)
      Hatom1 = find_protons(atom1, Coord)
      Hatom2 = find_protons(atom2, Coord)
      odl = []
      for mnum in range(1,21,1):
        Coor = eval('Coor' + str(mnum))
        dobs = getDistance(Hatom1, Hatom2, Coor)
        odl.append(dobs)
      observed.append(np.round(np.mean(odl),2))
      # if np.mean(odl) - upl < -2.0:
      #   print('{:} {:} mean {:3.2f} upl {:}'.format(atom1, atom2, np.mean(odl),upl))
      # if np.std(odl) > 1:
      #   print('{:} {:} std {:3.2f} upl {:}'.format(atom1, atom2, np.std(odl),upl))
  # fig=plt.figure(figsize=(4,4))
  # ax = fig.add_subplot(1,1,1)
  # ax.plot(cyana,cyana1,linewidth = 2, color = [0,0,0], label = None, zorder = 1)
  # ax.plot(cyana,cyana,linewidth = 2, color = [0,0,0], label = None, zorder = 1)
  # ax.scatter(cyana,observed, color = 'darkgreen',marker='o', s=30, label = None, clip_on=False, zorder = 2)
  # ax.set_xlabel('CYANA UPL')
  # ax.set_ylabel('Average Dist')
  # ax.set_ylim([1,8])
  # ax.set_ylim([1,8])
  # plt.show()

  DistancesDF.to_csv('post_cyana_analysis/' + pdb_name + '_distances_heavy_v3.csv')
  # print(DistancesDF)
  # print(DistancesDF.shape)
# print(DistancesDF.dropna(thresh=5).shape)
  return DistancesDF, FilteredDF
