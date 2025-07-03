import os
import sys
import pandas as pd
import numpy as np
import re
import math as m
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

#Atoms_dict = {'I':['CD1'], 'L':['CD1','CD2'], 'V':['CG1','CG2'], 'M':['CE'], 'A':['CB'], 'T':['CG2'], 'W':['NE1','HE1'], 'F':['CE1','CE2','HE1','HE2'], 'Y':['CE1','CE2','HE1','HE2'],'K':['CA','CB','CG','CD','CE','HA','HB3','HG3','HD3','HE2'],'R':['CA','CB','CG','CD','HA','HB3','HG3','HD3'],'H':['CE1','HE1','CD2','HD2']}
Atoms_dict = {'I':['CD1'], 'L':['CD1','CD2'], 'V':['CG1','CG2'], 'M':['CE'], 'A':['CB'], 'T':['CG2'],'F':['CE1','CE2'], 'Y':['CE1','CE2']}
AAA2A = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}


pdb1 = 
pdb2 = 

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

def find_atoms(inatom, PDBdict):
  atype = '{:}{:}'.format(inatom.split('-')[0][0], inatom.split('-')[-1])
  atoms = []
  if len(inatom.split('-')) == 2:
    resid, name = inatom.split('-')[0],inatom.split('-')[-1]
    for atomid in ['{:}-{:}'.format(resid,atom) for atom in Atoms_dict[resid[0]]]:
      if atomid in PDBdict.keys():
        atoms.append(atomid)
  if len(inatom.split('-')) == 3:
    resid, chain, name = inatom.split('-')[0],inatom.split('-')[1],inatom.split('-')[-1]
    for atomid in ['{:}-{:}-{:}'.format(resid,chain,atom) for atom in Atoms_dict[resid[0]]]:
      if atomid in PDBdict.keys():
        atoms.append(atomid)
  return atoms

for line in pdblines:

def get_coor(pdb, chain):
  Coor = {}
  Seq = []
  chainidx = 21
  if line[0:4] == "ATOM" and line[21] == ' ':
    chainidx = 72
    break 
  for line in pdb.readlines():
    if line == 'ENDMDL\n':
      break
    if line[0:4] == "ATOM" or line[0:4] == 'HETA':
      if line[17:20].strip() in AAA2A.keys() and if line[chainidx] == chain:
        if line[12:16].strip() in Atoms_dict[AAA2A[line[17:20].strip()]]:
          index = '{:}{:}-{:}'.format(AAA2A[line[17:20].strip()],line[22:26].strip(),line[12:16].strip())
          Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
          if '{:}{:}'.format(AAA2A[line[17:20].strip()],line[22:26].strip()) not in Seq:
            Seq.append('{:}{:}'.format(AAA2A[line[17:20].strip()],line[22:26].strip()))
  return Coor,Seq

Coor1, Seq1 = get_coor(pdb1, chain1)
Coor2, Seq2 = get_coor(pdb2, chain2)
dist1 = pd.DataFrame(index=Seq1, columns =Seq1)
dist2 = pd.DataFrame(index=Seq2, columns =Seq2)
ddm = pd.DataFrame(index=Seq1, columns =Seq1)

for res1 in Seq2:
  for res2 in Seq2:
    atom1 = find_atoms(res1,Coor2)
    atom2 = find_atoms(res2,Coor2)
    dist2.loc[res1,res2] = getDistance(atom1,atom2,Coor2)
for res1 in Seq1:
  m1a1 = find_atoms(res1,Coor1)
  m2a1 = find_atoms(res1,Coor2)
  if len(m1a1) !=0 and len(m2a1) != 0:
    for res2 in Seq1:
      if len(m1a2) !=0 and len(m2a2) != 0:
        m1a2 = find_atoms(res2,Coor1)
        m2a2 = find_atoms(res2,Coor2)
        dist1.loc[res1,res2] = getDistance(m1a1,m1a2,Coor1)
        dist2.loc[res1,res2] = getDistance(m2a1,m2a2,Coor2)
        ddm.loc[res1,res2] = round(dist1.loc[res1,res2] - dist2.loc[res1,res2])






