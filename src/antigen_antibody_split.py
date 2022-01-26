#!/Applications/anaconda3/envs/py27/bin/python
from pathlib import Path
import pymol
from pymol import cmd
import sys

# Pass the data path at the command line
dataPath = sys.argv[1]

PDBList = []
for PDB in dataPath.glob('PDBs/*_antigen_renumber.pdb'):
    filename = PDB.name
    filename_prefix = filename[:-4] # removes .pdb
    PDBList.append(filename_prefix)

print('{} PDBs being processed'.format(len(PDBList)))

for PDB in PDBList:
    print('  {} being procssed'.format(PDB))

    # Loads PDBs in PyMOL
    cmd.load('PDBs/'+pdb+'.pdb')

    # Selection strings to use in PyMOL
    # These could be updated
    antibody = "(m. " + pdb + " & c. H) + (m. " + pdb + " & c. L)"
    antigen = "(m. " + pdb + " & c. C) + (m. " + pdb + " & c. D) + (m. " + pdb + " & c. E) + (m. " + pdb + " & c. F) + (m. " + pdb + " & c. X) + (m. " + pdb + " & c. Y)"

    cmd.save('BSA/' + pdb + '/antibody.pdb', antibody)
    print('BSA/' + pdb + '/antibody.pdb')
    cmd.save('BSA/' + pdb + '/antigen.pdb', antigen)
    print('BSA/' + pdb + '/antigen.pdb')