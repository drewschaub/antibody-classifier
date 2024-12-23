import gzip
from pathlib import Path

def compress_pdb_files(input_dir, output_dir=None):
    """
    Compress all .pdb files in input_dir into .pdb.gz.
    By default, saves .pdb.gz files in the same directory unless output_dir is specified.
    """
    input_dir = Path(input_dir)
    if output_dir is None:
        output_dir = input_dir
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    for pdb_file in input_dir.glob("*.pdb"):
        output_file = output_dir / f"{pdb_file.stem}.pdb.gz"
        with open(pdb_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)
        print(f"Compressed {pdb_file.name} -> {output_file.name}")

if __name__ == "__main__":
    # Quick usage example
    compress_pdb_files("data/hiv1")
