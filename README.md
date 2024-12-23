Antibody Analysis Workflow
===================================

1. Download ```naccess``` from: http://www.bioinf.manchester.ac.uk/naccess/. Add ```naccess``` to your ```$PATH```



## Step 1 - Run naccess and analyze output
1. Prepare PDBs
	Suggestions:
	* Unify chain labeling of antibody and antigen chains (e.g. H/L for heavy and light chains, respectively)
	* One antibody per antigen
	* Save structure as two PDB files:
		* antigen.pdb containing only antigen
		* antibody.pdb containing only antibody

		* Automate splitting PDB of complex structure into antigen.pdb and antibody.pdb using:

				antigen_antibody_split.py

2. Run naccess:

		get_asa.py

	to obtain accessible surface area of antibody, antigen, and antibody-antigen complex.

3. Parse naccess output

	Run all three python scripts in the naccess folder. For class ID, only calculate_epitope_residues.py is necessary. The 2 paratope pythons will be for identifying conserved ab-antigen interactions for sub-epitope classes.

	*  Calculate the epitope residues
		
			calculate_epitope_residues.py
	
		Output: `filename_epitope_none_zero_bsa.out`

	* Calculate the paratope residues
	
			calculate_paratope_residues.py
	
		Output: `filename_paratope_none_zero_bsa.out`
	* Parse the paratope bsa file
	
			parse_paratope_bsa.py
      

# Curating Datasets

Store your pdb files in `pdb.gz` format in the `data` folder. You can separate projects. To match with the manuscript this has folders for `hiv1`, `influenza` and `sarscov2`. To gzip on MacOS simply copy your `pdb` files to your folder and run `gzip *.pdb`. It will automatically compress to `pdb.gz` and remove the orignal `pdb` files. These are usually ~20% of the original size. 

All antibodies should use the same naming convention. We recommend `IMGT` as that is what was used for the manuscript. Truncate antibody chains after FR4, or after res 128 in IMGT numbering.