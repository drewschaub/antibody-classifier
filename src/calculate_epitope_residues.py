#!/usr/bin/python
#
# Original code from Jing Zhou (Dec 15, 2019)
# Updated by Yiran Wang (Mar 29, 2020)
# Updated by Andrew Schaub (Jan 26, 2022) 
# Calculate epitope buried surface 
#
# Action Items: 
#  
#   1. I'm not comfortable with the os system calls
#   2. Don't need to use file, we can use pathlib, then we can ignore
#      all the opens and closes as pathlib will take care of that for us
from pathlib import Path
import os

target = 'influenza'
targetDataPath = Path('../../../data/{}/BSA/tmp'.format{target})
naccessPath = 'naccess'

for file in targetDataPath.glob('*'):
	print('Processing: {}'.format(file))

	# systems calls are dangerous I need to update this to something secure
	os.system("cd {0} && {}/combine_rsa.pl antigen.rsa antigen_antibody.rsa > epitope.out".format(file, naccessPath))
	os.system("cd {0} && {}/get_glycan_per.pl epitope.out > epitope_bsa.out".format(file, naccessPath))

	# This needs updated to use paths. 
	#Analyze results
	inputFile = open('{}/epitope.out'.format(file), 'r')
	outputFile = open('{}/{}_epitope_none_zero_bsa.out'.format(file, file), 'w')

	for line in inputFile:
		surface = float(line.split()[-1])
		if surface != 0:
			outputFile.write(line) # can use write_path here with pathlib
	outputFile.close()