import numpy as np
from pymol import cmd
from .residue_ops import flip_sidechain_if_ambiguous
import scipy.cluster.hierarchy as sch

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

def compute_rmsd_single(mob, targ, chain_type, epitope_selection, antigen_chain="W"):
    """Compute RMSD between two structures for a specific chain"""
    # Align the structures
    print(f"Aligning {mob}_{chain_type} to {targ}_{chain_type}")
    cmd.align(f"{mob}_{chain_type}", f"{targ}_{chain_type}")
    
    # Create selections
    mobile = f"{mob} and chain {antigen_chain} and i. {epitope_selection}"
    target = f"{targ} and chain {antigen_chain} and i. {epitope_selection}"
    print(f"Mobile: {mobile}, Target: {target}")

    # Get initial RMSD
    rmsd = cmd.rms_cur(mobile, target, matchmaker=-1)
    print(f"Initial RMSD: {rmsd}")
    
    # Check for ambiguous residues
    atoms = cmd.get_model(mobile)
    try:
        if atoms.atom[0].resn in AMBIGUOUS_RESIDUES:
            flipped = flip_sidechain_if_ambiguous(f"chain {antigen_chain} and i. {epitope_selection}")
            rmsd_flipped = cmd.rms_cur(mobile, flipped)
            rmsd = min(rmsd, rmsd_flipped)
    except IndexError:
        print(f"Index Error: {mob} and chain {antigen_chain}")
        
    return rmsd

def compute_rmsd_matrix(structure_list, chain_mode='both', epitope_selection='all'):
    """
    Compute RMSD matrix between structures.
    chain_mode: 'both', 'H', or 'L'
    """
    n = len(structure_list)
    rmsd_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Validate chain_mode
            if chain_mode not in ['both', 'H', 'L']:
                raise ValueError(f"Invalid chain_mode: {chain_mode}")
                
            # Compute RMSD based on chain mode
            chain_type = 'HL' if chain_mode == 'both' else chain_mode
            rmsd = compute_rmsd_single(
                structure_list[i],
                structure_list[j], 
                chain_type,
                epitope_selection
            )
            
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
