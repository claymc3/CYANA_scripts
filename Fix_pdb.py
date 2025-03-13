import os
import sys
import re
import pandas as pd
import numpy as np

inpdb = sys.argv[1]
correction = float(sys.argv[2])


if len(sys.argv)==1:
  print('''

Usage: 
    fixpdb [PDB] [numbering correction]

Required Input:

    PDB         assuming your running script in CNS/refinePDB/ directory 
                name is the name you supplied to the cnsprep and cnspost
                scripts 

    numbering correction    how much should be added to correct the numbering

OutPut:
    PDB_fix.tbl

''')

newpdb = open(inpdb.replace('.pdb','_fixed.pdb'),'w')
# for line in open(inpdb).readlines():
#   if re.match(r"^ATOM", line):
#       resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)
#       newline = "{:} {:} {:<4}{:}".format(line[0:20],line[72],resid,line[27:])
#       newpdb.write(newline)
#   if not re.match(r"^ATOM", line):
#       newpdb.write(line)
# newpdb.close()
for line in open(inpdb).readlines():
  if re.match(r"PTR", line[17:20].strip()):
    resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)
    line1 = line.replace('ATOM  ','HETATM').replace('1HE ',' HE1').replace('2HE ',' HE2').replace('1HD ',' HD1').replace('2HD ',' HD2')
    resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)
    newline = "{:} {:} {:<4}{:}".format(line1[0:20],line[72],resid,line1[27:])
    # newline = "{:} {:}{:}".format(line1[0:20],line[72],line1[22:])
    newpdb.write(newline)
  if re.match(r"^ATOM", line) and not re.match(r"PTR", line[17:20].strip()):
    newline = "{:} {:}{:}".format(line[0:20],line[72],line[22:])
    newpdb.write(newline)
  if not re.match(r"^ATOM", line):
    newpdb.write(line)
newpdb.close()

  # if re.match(r"TPO", line[17:20].strip()):
  #   print(line)
  #   resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)
  #   line1 = line.replace('ATOM  ','HETATM').replace('1HG2','HG21').replace('2HG2','HG22').replace('3HG2','HG23')
  #   newline = "{:} {:} {:<4}{:}".format(line1[0:20],line[72],resid,line1[27:])
  #   newpdb.write(newline)

  # if re.match(r"SEP", line[17:20].strip()):
  #   print(line)
  #   resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)
  #   line1 = line.replace('ATOM  ','HETATM')
  #   newline = "{:} {:} {:<4}{:}".format(line1[0:20],line[72],resid,line1[27:])
  #   newpdb.write(newline)