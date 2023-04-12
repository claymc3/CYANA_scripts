"""
Originally written by Mary Clay
St. Jude Children's Research Hospital
April 04, 2023, MC

Changelog:

"""

import os
import sys
import re

if len(sys.argv) < 3:
	print('''

		Alias this script in the appropriate shell script.

		Usage:

		findpeaks residue_1 residue_2

		residue_1
			residue idenity (L481)
			sequence index (481)
			resoance name (L481-QD1)

		residue_2
			residue idenity (M518)
			sequence index (518)
			resoance name (M518-H)


		Examples (assuming findpeaks is an alias to this script):

		findpeaks 481 518
			- Returns matching entries containing 481 and 518
			- order is ignored
			- Empty if no matches exist

		findpeaks L481 M518
			- Returns matching entries containing L481 and M518
			- order is ignored
			- Empty if no matches exist
		
		findpeaks L481-QD1 M518
			- Returns matching entries containing L481-QD1 and any M518
			- order is ignored
			- Empty if no matches exist
		
		uplsearch final.upl L481-QD1 M518-QE
			- Returns matching entries containing L481-QD1 and any M518-QE
			- order is ignored
			- Empty if no matches

		''')
	exit()



# First res/atom specs
ra1 = sys.argv[1]

# Second res/atom specs
ra2 = sys.argv[2]
regex1 = r".*"+re.escape(ra1)+r".*"+re.escape(ra2)+r".*"
regex2 = r".*"+re.escape(ra2)+r".*"+re.escape(ra1)+r".*"

for line in open('noa_analysis/Assignment_Summary.txt').readlines():
	if line.strip():
		selem = line.replace('#','').split()[0]
		if re.findall(regex1, selem ) or re.findall(regex2, selem ):
			print(line[:-1])

