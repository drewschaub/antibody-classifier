#!/usr/bin/python
#
# I need to update these old scripts. I don't use os or glob anymore
# as of Python 3.4 I prefer pathlib
# Ref: https://docs.python.org/3/library/pathlib.html
import sys, os, glob

naccessPath = 'naccess2.1.1/naccess'

#folder = ['A12V163-b.01_trimer']
#for i in folder:
for file in glob.glob('influenza/*'):
	print(file)
    
	cmd1 = 'cat ' + file +'/antigen.pdb '  + file + '/antibody.pdb | grep -v "END" > ' + file +'/antigen_antibody.pdb'
	print(cmd1)

	os.system(cmd1)
	os.system("cd {0} && {} -h antibody.pdb".format(file, naccessPath))
	os.system("cd {0} && {} -h antigen.pdb".format(file, naccessPath))
	os.system("cd {0} && {} -h antigen_antibody.pdb".format(file, naccessPath))
