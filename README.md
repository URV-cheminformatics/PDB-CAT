# PDB-CAT
<img src="image_documentation/PDB-CAT.jpg" width="300">

PDB files contain structural information about proteins and other biomolecules, and are widely used in Drug Discovery. However, sorting through large numbers of PDB files to find the desired structures can be a time-consuming task. PDB-CAT is a Jupyter Notebook that aims to simplify this process by automatically categorizing the structures based on the type of interaction between atoms in the protein and the ligand, and checking for any mutations in the sequence. PDB-CAT is a program that classifies a group of protein structures into three categories: covalent bonded, non-covalent bonded, and no-bonded. Before classification, the program verifies if there are any mutations in the protein sequence by searching for the amino acid sequence located in the SEQRES lines of the PDB file and comparing it to a reference sequence. The program outputs a file with structures that meet the specified criteria, which can then be used for further analysis, such as virtual screening or molecular docking. The program is easy to use and can be customized to fit the user's specific needs.

## REQUIREMENTS
This notebook uses Python 3.10.9
This program requires the following packages:
- Biopython (1.81)
- os
- re
- csv
- shutil
- pdbecif

You can install these packages using pip:
    $ pip install <library>

*You do not need extensive knowledge of these libraries.*

## EXPLANATION
<img src="image_documentation/PDB-CAT_poster.png" width="300">

