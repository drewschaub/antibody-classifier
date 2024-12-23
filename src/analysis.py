import numpy as np
from pymol import cmd
from .residue_ops import switch_atom_name
from concurrent.futures import ThreadPoolExecutor, as_completed
import scipy.cluster.hierarchy as sch
from pathlib import Path
import json
import random
import pandas as pd

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

    unique_mob_name = generate_unique_name(mob)
    unique_targ_name = generate_unique_name(targ)

    cmd.create(unique_mob_name, mob)
    cmd.create(unique_targ_name, targ)

    # generate selection strings
    chains = {'HL': 'chain H+L',
              'H': 'chain H',
              'L': 'chain L'}
    
    mob_selection = f"{unique_mob_name} and {chains[chain_type]} and backbone"
    targ_selection = f"{unique_targ_name} and {chains[chain_type]} and backbone"

    # align structures
    try:
        # print(f"align {mob_selection}, {targ_selection}")
        cmd.align(mob_selection, targ_selection)
    except:
        # output a PSE file to look at the error
        cmd.save(f"{mob}_{chain_type}_error.pse")
        # print(f"Error aligning {mob}_{chain_type} and {targ}_{chain_type}")

    # Calculate RMSD
    mob_epitope_sele = f"{unique_mob_name} and {epitope_selection}"
    targ_epitope_sele = f"{unique_targ_name} and {epitope_selection}"
    # print(mob_epitope_sele, targ_epitope_sele)
    rmsd = cmd.rms_cur(mob_epitope_sele, targ_epitope_sele)

    # Create a copy of the mobile object
    flipped_mob = generate_unique_name(f"{mob}_flipped")

    cmd.create(flipped_mob, unique_mob_name)
    
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
                    # print(rmsd)
                else:
                    # switch back to original atom names
                    switch_atom_name(residue_selection, flip[1], flip[0])
   
    # delete the temporary objects
    cmd.delete(flipped_mob)
    cmd.delete(unique_mob_name)
    cmd.delete(unique_targ_name)

    return rmsd

def compute_rmsd_matrix(data_dir, structure_list, chain_mode='HL', epitope_selection='all'):
    """
    Compute RMSD matrix between structures.
    chain_mode: 'HL', 'H', or 'L'
    """
    n = len(structure_list)
    rmsd_matrix = np.zeros((n, n))
    new_antibodies = structure_list

    # read 'epitope_dict.csv' and convert to epitope selection string
    epitope_dict_path = Path(data_dir, 'epitope_dict.csv')
    epitope_dict = {}

    if epitope_dict_path.exists():
        epitope_df = pd.read_csv(epitope_dict_path)
        epitope_dict = dict(zip(epitope_df['structure_stem'], epitope_df['epitope_dict']))
            
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
    
    if chain_mode not in ['HL', 'H', 'L']:
        raise ValueError(f"Invalid chain_mode: {chain_mode}")

    def compute_pairwise_rmsd(i, j):
        chain_type = chain_mode
        if structure_list[i] in new_antibodies or structure_list[j] in new_antibodies:
            epitope_dict_i = parse_epitope_dict(epitope_dict.get(structure_list[i], '{}'))
            epitope_dict_j = parse_epitope_dict(epitope_dict.get(structure_list[j], '{}'))
            # print(epitope_dict_i)
            # print(epitope_dict_j)
            intersection_dict = find_epitope_intersection(epitope_dict_i, epitope_dict_j)
            epitope_selection = f"({generate_epitope_selection(intersection_dict)})"
            # print(epitope_selection)
            return (i, j, compute_rmsd_single(structure_list[i], structure_list[j], chain_type, epitope_selection))
        else:
            old_i = old_structure_list.index(structure_list[i])
            old_j = old_structure_list.index(structure_list[j])
            return (i, j, old_rmsd_matrix[old_i, old_j])

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(compute_pairwise_rmsd, i, j) for i in range(n) for j in range(i+1, n)]
        for future in as_completed(futures):
            i, j, rmsd = future.result()
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
            
    return rmsd_matrix

def generate_unique_name(base_name):
    while True:
        random_int = random.randint(0, 10000)
        unique_name = f"{base_name}_{str(random_int)}"
        if unique_name not in cmd.get_names("objects"):
            return unique_name

def hierarchical_clustering(matrix, method='average'):
    """
    Perform hierarchical clustering on the RMSD matrix.
    Return the linkage matrix, etc.
    """
    Y = sch.linkage(matrix, method=method)
    return Y

def parse_epitope_dict(epitope_dict_str):
    epitope_dict = json.loads(epitope_dict_str.replace("'", "\""))
    return epitope_dict

def find_epitope_intersection(epitope_dict_i, epitope_dict_j):
    intersection_dict = {}
    for chain in epitope_dict_i.keys() & epitope_dict_j.keys():
        residues_i = set(epitope_dict_i[chain])
        residues_j = set(epitope_dict_j[chain])
        intersection_residues = residues_i & residues_j
        if intersection_residues:
            intersection_dict[chain] = list(intersection_residues)
    return intersection_dict

def generate_epitope_selection(epitope_dict):
    selections = []
    for chain, residues in epitope_dict.items():
        resi_str = "+".join(map(str, residues))
        selections.append(f"(chain {chain} and resi {resi_str})")
    return "+".join(selections)