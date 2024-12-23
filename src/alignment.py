from pymol import cmd
from .residue_ops import flip_sidechain_if_ambiguous

def align_structures(ref_selection, mob_selection):
    """
    Perform a PyMOL 'align' command between two selections.
    Return a tuple (rmsd_value, alignment_info).
    """
    try:
        alignment_data = cmd.align(ref_selection, mob_selection)
        # alignment_data is usually a tuple from PyMOL: (RMSD, ...), etc.
        rmsd_value = alignment_data[0]
        return rmsd_value, alignment_data
    except Exception as e:
        print(f"Error performing alignment: {e}")
        return None, None

def measure_rms_cur(ref_selection, mob_selection):
    """
    After an alignment, you can do cmd.rms_cur if you want a separate measure. 
    """
    try:
        return cmd.rms_cur(ref_selection, mob_selection)
    except Exception as e:
        print(f"Error computing RMS: {e}")
        return None

def align_with_sidechain_flip(ref_selection, mob_selection):
    """
    Illustrates how you might incorporate flipping logic if needed.
    """
    base_rms = measure_rms_cur(ref_selection, mob_selection)
    # Example of flipping logic on the mobile selection:
    flipped_sel = flip_sidechain_if_ambiguous(mob_selection)
    flipped_rms = measure_rms_cur(ref_selection, flipped_sel)
    return min(base_rms, flipped_rms)
