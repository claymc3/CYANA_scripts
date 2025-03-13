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
    fixpdb [name] [numbering correction]

Required Input:

    name        assuming your running script in CNS/refinePDB/ directory 
                name is the name you supplied to the cnsprep and cnspost
                scripts 

    numbering correction    how much should be added to correct the numbering

OutPut:
    name_noe_rn.tbl
    name_dihed_rn.tbl
    name_hbond_rn.tbl
''')

newpdb = open(inpdb.replace('.pdb','_fixed.pdb'),'w')
for line in open(inpdb).readlines():
	if re.match(r"^ATOM", line):
		resid = '{:<2.0f}'.format(float(line[22:26].strip()) + correction)

		newline = "{:} {:} {:<4}{:}".format(line[0:20],line[72],resid,line[27:])
		newpdb.write(newline)
	if not re.match(r"^ATOM", line):
		newpdb.write(line)
newpdb.close()