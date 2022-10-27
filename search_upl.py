"""
Originally written by Mary Clay
St. Jude Children's Research Hospital

Modified by Mukundan Ragavan
St. Jude Children's Research Hospital



Changelog:

Sep 28, 2022, MR
 - add additional examples

Oct 03, 2022, MR
 - improve regex; add new conditionals and add aa_list

"""

import os
import sys
import re

if len(sys.argv) < 3:
	print('''

		Alias this script in the appropriate shell script.

		Usage:

		uplsearch <upl file(s)> <sequence position of amino acid> [optional: <atom code> or <sequence position>]

		upl file(s): One or more upl files can be speicified - comma separated

		sequence position used in assignments

		Optional argument:
		1- to 3- letter atom or residue identifier (e.g. H, QD1, QE, 726, etc.) or sequence position. 
		atom or residue identifier is case insensitive.


		Examples (assuming uplsearch is an alias to this script):

		uplsearch final.upl 688 qg2
			- Returns matching entries containing 688 QG2
			- Empty if no matches exist

		uplsearch final.upl 688 val
			- Returns matching entries containing 688 QG2
			- Empty if no matches exist
		
		uplsearch final.upl 688
			- Returns matching entries containing 688
			- Empty if no matches exist
		
		uplsearch final.upl 688 726
			- Returns matching entries containing 688 and 726
			- order is ignored
			- Empty if no matches


		''')
	exit()


aa_list = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']

# list of filenames separated by comma
sys.argv[1]

# residue number
res_num = sys.argv[2]

# Atom identifier - 1- to 3- letter
# case insensitive
if len(sys.argv) == 4:
	opt_arg = sys.argv[3]

for fnam in sys.argv[1].split(','):
	print('\n*** upl file : ' + fnam + '\n')
	for line in open(fnam).readlines():
		if len(sys.argv) == 4:
			if opt_arg.isdigit():
				regex1 = r".*"+re.escape(res_num)+r"\s*[A-Z]{3}.*"+re.escape(opt_arg)+r".*"
				regex2 = r".*"+re.escape(opt_arg)+r"\s*[A-Z]{3}.*"+re.escape(res_num)+r".*"

				if re.findall(regex1, line) or re.findall(regex2, line):
					print(line.strip())

			else:
				if opt_arg.upper() in aa_list:
					regex = r".*"+re.escape(res_num)+r"\s"+re.escape(opt_arg.upper())
				else:
					regex = r".*"+re.escape(res_num)+r"\s*[A-Z]{3}\s{2}"+re.escape(opt_arg.upper() + ' ')+r".*"

				if re.findall(regex, line):
					print(line.strip())

		else:
			regex = r".*"+re.escape(res_num)+r"\s*[A-Z]{3}.*"

			if re.findall(regex, line):
				print(line.strip())
