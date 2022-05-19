import os
import sys

sys.argv[1]
sys.argv[2]
sys.argv[3]
for line in open(sys.argv[1]).readlines():
	if sys.argv[2] + ' ' + sys.argv[3] in line:
		print(line.strip())