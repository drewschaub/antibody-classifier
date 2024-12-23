import argparse
import sys
from pathlib import Path
import pandas as pd
from pymol import cmd
import json
import numpy as np
from .structure_io import init_pymol_headless, load_structure, create_chain_selection
# from .data_prep import compress_pdb_files
# from .alignment import align_structures, measure_rms_cur
from .analysis import compute_rmsd_matrix, hierarchical_clustering
# from .visualization import plot_rmsd_heatmap, plot_dendrogram

def main():
    parser = argparse.ArgumentParser(description="Antibody RMSD alignment pipeline.")
    parser.add_argument("--data-dir", default="data/hiv1",
                        help="Path to input data with PDB files")
    args = parser.parse_args()

    init_pymol_headless()
    
    # Example of loading structures
    data_dir = Path(args.data_dir)
    pdb_files = list(data_dir.glob("*.pdb.gz")) + list(data_dir.glob("*.pdb"))
    # sort pdb files by name
    pdb_files.sort()
    
    structure_names = []
    for pdb_file in pdb_files:
        obj_name = pdb_file.stem  # strip .pdb or .pdb.gz
        load_structure(pdb_file, obj_name)
        structure_names.append(obj_name)
    
    # Create chain selections or do other manipulations
    for name in structure_names:
        create_chain_selection(name, 'H', f"{name}_H", backbone_only=True)
        create_chain_selection(name, 'L', f"{name}_L", backbone_only=True)
        create_chain_selection(name, 'HL',f"{name}_HL", backbone_only=True)
    
    # save pse file
    cmd.save("output.pse")
    
    # we'll use a temporary epitope selection until i bring in the freesasa information
    epitope_selection = 'c. B and i. 512-519'

    # compute the RMSD matrices for H chains
    rmsd_matrix_H = compute_rmsd_matrix(data_dir, structure_names, chain_mode='H', epitope_selection=epitope_selection)
    cmd.save("rmsd_matrix_H.pse")

    # save as a csv where row names and column names are the structure names
    rmsd_df_H = pd.DataFrame(rmsd_matrix_H, index=structure_names, columns=structure_names)
    rmsd_df_H.to_csv(Path(data_dir, "rmsd_matrix_H.csv"))

    # save metadata as json
    metadata = {"row_labels": structure_names, "col_labels": structure_names}
    metadata_path = Path(data_dir, "rmsd_matrix_H.json")
    with open(metadata_path, "w") as f:
        json.dump(metadata, f)

    # save numpy array
    rmsd_matrix_H_path = Path(data_dir, "rmsd_matrix_H.npy")
    np.save(rmsd_matrix_H_path, rmsd_matrix_H)




    # Example RMSD matrix computation
    # (You'd need to adapt compute_rmsd_matrix to do the actual alignment calls)
    # rmsd_matrix = compute_rmsd_matrix(structure_names, chain_mode='H_noH')
    
    # Visualize
    # plot_rmsd_heatmap(rmsd_matrix, structure_names, title="Heavy-Chain RMSD Heatmap")
    # ...
    
    # Save or show
    # ...
    return 0

if __name__ == "__main__":
    sys.exit(main())
