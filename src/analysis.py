import numpy as np
from pymol import cmd
from .residue_ops import switch_atom_name
import scipy.cluster.hierarchy as sch
from pathlib import Path
import json

# Constants
AMBIGUOUS_RESIDUES = ['HIS',
                      'ASP',
                      'ARG',
                      'PHE',
                      'GLN',
                      'GLU',
                      'LEU',
                      'ASN',
                      'TYR',
                      'VAL'
                      ]

def compute_rmsd_single(mob, targ, chain_type, epitope_selection):
    """Compute RMSD between two structures for a specific chain"""

    rmsd = 9999.99

    # align structures using antibody
    try:
        print(f"align {mob}_{chain_type} and backbone, {targ}_{chain_type} and backbone")
        cmd.align(f"{mob}_{chain_type} and backbone", f"{targ}_{chain_type} and backbone")
    except:
        # output a PSE file to look at the error
        cmd.save(f"{mob}_{chain_type}_error.pse")
        # print(f"Error aligning {mob}_{chain_type} and {targ}_{chain_type}")

    
    # Create epitope selections
    mob_epitope_sele = f"{mob} and {epitope_selection}"
    targ_epitope_sele = f"{targ} and {epitope_selection}"
    print(f"sele mobile, {mob_epitope_sele}")
    print(f"sele target, {targ_epitope_sele}")

    # Create a copy of the mobile object
    flipped_mob = f"{mob}_flipped"
    cmd.create(flipped_mob, mob)
    
    flips = {'ARG': [('NH1', 'NH2')],
             'HIS': [('ND1', 'CD2'), ('ND1', 'CD2')],
             'ASP': [('OD1', 'OD2')],
             'PHE': [('CD1', 'CD2'), ('CE1', 'CE2')],
             'GLN': [('OE1', 'NE2')],
             'GLU': [('OE1', 'OE2')],
             'LEU': [('CD1', 'CD2')],
             'ASN': [('OD1', 'ND2')],
             'TYR': [('CD1', 'CD2'),('CE1', 'CE2')],
             'VAL': [('CG1', 'CG2')]}

    # Iterate through the epitope selection alpha carbons
    atoms = cmd.get_model(mob_epitope_sele)
    for atom in atoms.atom:
        if atom.resn in AMBIGUOUS_RESIDUES:
            res_flips = flips[atom.resn]
            for flip in res_flips:
                # select the correct chain and residue to flip on mobile
                residue_selection = f"{flipped_mob} and chain {atom.chain} and resi {atom.resi}"
                switch_atom_name(residue_selection, flip[0], flip[1])

                # calculate rmsd of flipped res obj
                flipped_rmsd = cmd.rms_cur(f"{flipped_mob} and {epitope_selection}", f"{targ} and {epitope_selection}")
                
                # if the flipped rmsd is less than the original rmsd, keep the flipped residue
                if flipped_rmsd < rmsd:
                    rmsd = flipped_rmsd
                    print(rmsd)
                else:
                    # switch back to original atom names
                    switch_atom_name(residue_selection, flip[1], flip[0])
   
    return rmsd

def compute_rmsd_matrix(data_dir, structure_list, chain_mode='both', epitope_selection='all'):
    """
    Compute RMSD matrix between structures.
    chain_mode: 'both', 'H', or 'L'
    """
    n = len(structure_list)
    rmsd_matrix = np.zeros((n, n))

    new_antibodies = structure_list

    # check data_dir and see if 
    old_rmsd_matrix_np_path = Path(data_dir, f"rmsd_matrix_{chain_mode}.npy")
    old_rmsd_matrix_json_path = Path(data_dir, f"rmsd_matrix_{chain_mode}.json")

    if old_rmsd_matrix_np_path.exists() and old_rmsd_matrix_json_path.exists():
        print(f"Loading old RMSD matrix from {old_rmsd_matrix_np_path}")
        old_rmsd_matrix = np.load(old_rmsd_matrix_np_path)
        with open(old_rmsd_matrix_json_path, "r") as f:
            metadata = json.load(f)
            old_structure_list = metadata['row_labels']
        
            if old_structure_list == structure_list:
                print("No new antibodies added, using old RMSD matrix")
                return old_rmsd_matrix
            else:
                # get list of new antibodies
                new_antibodies = list(set(structure_list) - set(old_structure_list))
        
    for i in range(n):
        for j in range(i+1, n):
            print(i, j)
            # Validate chain_mode
            if chain_mode not in ['both', 'H', 'L']:
                raise ValueError(f"Invalid chain_mode: {chain_mode}")
                
            # Compute RMSD based on chain mode
            chain_type = 'HL' if chain_mode == 'both' else chain_mode

            # if structure_list[i] or structure_list[j] in new_antibodies, compute rmsd, else use old rmsd
            if structure_list[i] in new_antibodies or structure_list[j] in new_antibodies:
                rmsd = compute_rmsd_single(
                    structure_list[i],
                    structure_list[j], 
                    chain_type,
                    epitope_selection
                )
            else:
                print('found old rmsd')
                # get index values of structure_list[i] and structure_list[j] in old_structure_list
                old_i = old_structure_list.index(structure_list[i])
                old_j = old_structure_list.index(structure_list[j])
                rmsd = old_rmsd_matrix[old_i, old_j]
            
            # Fill both sides of symmetric matrix
            rmsd_matrix[i,j] = rmsd
            rmsd_matrix[j,i] = rmsd
            
    return rmsd_matrix

def hierarchical_clustering(matrix, method='average'):
    """
    Perform hierarchical clustering on the RMSD matrix.
    Return the linkage matrix, etc.
    """
    Y = sch.linkage(matrix, method=method)
    return Y
