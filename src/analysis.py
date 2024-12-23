import numpy as np
import scipy.cluster.hierarchy as sch

def compute_rmsd_matrix(structure_list, chain_mode='both'):
    """
    Example function that loops over structures,
    does alignment, and populates a 2D RMSD array.
    chain_mode could be 'both', 'H', 'L', etc.
    """
    n = len(structure_list)
    rmsd_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Insert code to do alignment in PyMOL, measure RMSD
            # e.g. align_structures(f"{structure_list[i]}_{chain_mode}", 
            #                       f"{structure_list[j]}_{chain_mode}")
            # ...
            # rmsd_matrix[i,j] = ...
            # rmsd_matrix[j,i] = ...
            pass
    
    return rmsd_matrix

def hierarchical_clustering(matrix, method='average'):
    """
    Perform hierarchical clustering on the RMSD matrix.
    Return the linkage matrix, etc.
    """
    Y = sch.linkage(matrix, method=method)
    return Y
