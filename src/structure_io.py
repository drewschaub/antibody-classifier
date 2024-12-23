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
        object_name = structure_path.stem

    sanitized_name = sanitize_name(object_name)
    
    try:
        cmd.load(structure_path.as_posix(), sanitized_name)
        print(f"Loaded structure: {sanitized_name} from {structure_path}")
        return sanitized_name # Return santiized name for further use
    except Exception as e:
        print(f"Error loading {structure_path}: {e}")
        return None

def create_chain_selection(object_name, chain_id, sel_name=None, backbone_only=True):
    """
    Create a PyMOL selection for a given chain.
    Sanitize names for PyMOL compatibility.
    """
    if sel_name is None:
        sel_name = f"{object_name}_{chain_id}_sel"
    
    sel_name = sanitize_name(sel_name)
    object_name = sanitize_name(object_name)
    
    selection = f"(m. {object_name} & c. {chain_id})"

    if backbone_only:
        selection += " & n. CA+C+N+O"
    
    try:
        cmd.select(sel_name, selection)
        print(f"Selection {sel_name} created: {selection}")
        return sel_name
    except Exception as e:
        print(f"Error creating selection {sel_name}: {e}")
        return None

def sanitize_name(name):
    """
    Sanitize a name for use in PyMOL (remove or replace invalid characters).
    """
    sanitized = name.replace('[', '').replace(']', '')
    return sanitized
