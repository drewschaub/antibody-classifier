# residue_ops.py

"""
This code is partially adapted from:
https://pymolwiki.org/index.php/RmsdByResidue
"""
from pymol import cmd

def switch_atom_name(residue_selection, atom_name1, atom_name2):
    cmd.alter(f"{residue_selection} and name {atom_name1}", 'name="tempAtom"')
    cmd.alter(f"{residue_selection} and name {atom_name2}", f'name="{atom_name1}"')
    cmd.alter(f"{residue_selection} and name tempAtom", f'name="{atom_name2}"')
    cmd.sort()

def flip_sidechain_if_ambiguous(residue_selection):
    """
    Flips atoms in ambiguous residues, returning the new object name.
    """
    flips = {
        "ASP": ("OD1", "OD2"),
        "ASN": ("OD1", "ND2"),
        "GLU": ("OE1", "OE2"),
        "GLN": ("OE1", "NE2"),
        "ARG": ("NH1", "NH2"),
        "LEU": ("CD1", "CD2"),
        "HIS": ("ND1", "NE2"),
        "TYR": ("CD1", "CD2"),
        "PHE": ("CD1", "CD2"),
        # "VAL" might need a different logic or none
    }

    flipped_name = "flippedRes"
    cmd.create(flipped_name, residue_selection)

    # Iterate atoms in the new object
    model = cmd.get_model(flipped_name)
    for atom in model.atom:
        if atom.resn in flips:
            switch_atom_name(
                f"{flipped_name} and resi {atom.resi} and resn {atom.resn}",
                flips[atom.resn][0],
                flips[atom.resn][1]
            )
    return flipped_name