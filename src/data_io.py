# src/antibody_classifier/data_io.py
import gzip
import os
from Bio.PDB import PDBParser

def load_structures_from_directory(directory_path):
    """
    Loads all .pdb or .pdb.gz files from a directory using Biopython.
    
    :param directory_path: Path to a folder containing PDB or PDB.GZ files.
    :return: Dict[str, Bio.PDB.Structure.Structure]
    """
    parser = PDBParser(QUIET=True)
    structures = {}
    
    for filename in os.listdir(directory_path):
        filepath = os.path.join(directory_path, filename)
        if filename.lower().endswith(".pdb.gz"):
            structure_id = filename[:-7]  # remove .pdb.gz
            with gzip.open(filepath, "rt") as handle:
                # parse structure
                structure = parser.get_structure(structure_id, handle)
                structures[structure_id] = structure
        elif filename.lower().endswith(".pdb"):
            structure_id = filename[:-4]  # remove .pdb
            structure = parser.get_structure(structure_id, filepath)
            structures[structure_id] = structure
    
    return structures
