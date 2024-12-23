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

    rmsd_matrices = []
    for aligntype in ['H', 'L', 'HL']:
        print(f"Computing RMSD matrix for {aligntype} chains")
        rmsd_matrix = compute_rmsd_matrix(data_dir, structure_names, chain_mode=aligntype, epitope_selection=epitope_selection)
        
        # save as a csv where row names and column names are the structure names
        rmsd_df = pd.DataFrame(rmsd_matrix, index=structure_names, columns=structure_names)
        rmsd_df.to_csv(Path(data_dir, f"rmsd_matrix_{aligntype}.csv"))

        # save metadata as json
        metadata = {"row_labels": structure_names, "col_labels": structure_names}
        with open(Path(data_dir, f"rmsd_matrix_{aligntype}.json"), "w") as f:
            json.dump(metadata, f)

        # save numpy array
        np.save(Path(data_dir, f"rmsd_matrix_{aligntype}.npy"), rmsd_matrix)

        rmsd_matrices.append(rmsd_matrix)



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
