import pymol
from pymol import cmd
from pathlib import Path

def init_pymol_headless():
    """
    Initialize PyMOL in headless mode if needed.
    """
    import __main__
    __main__.pymol_argv = ['pymol','-cqi']  # or other flags you prefer
    pymol.finish_launching()

def load_structure(structure_path, object_name=None):
    """
    Load a .pdb or .pdb.gz file into PyMOL.
    If object_name is not given, use file stem as the object name.
    """
    structure_path = Path(structure_path)
    if object_name is None:
        # object_name should not contain '.pdb.gz' or '.pdb'
        object_name = structure_path.stem
        object_name = object_name.replace('.pdb', '').replace('.pdb.gz', '')
    
    try:
        cmd.load(structure_path.as_posix(), object_name)
        print(f"Loaded structure: {object_name} from {structure_path}")
    except Exception as e:
        print(f"Error loading {structure_path}: {e}")

def create_chain_selection(structure_name, chain_id, sel_name=None, backbone_only=True):
    """
    Create a PyMOL selection for a given chain.
    e.g. create_chain_selection('myStruct', 'H', 'myStruct_H')
    """
    if sel_name is None:
        sel_name = f"{structure_name}_{chain_id}_sel"
    selection = f"(m. {structure_name} & c. {chain_id})"
    cmd.select(sel_name, selection)
    print(f"Selection {sel_name} defined as: {selection}")
    return sel_name
