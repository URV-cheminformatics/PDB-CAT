# PDB-CAT
![Image URL](image_documentation/PDB-CAT.png | width = 800)

PDB files contain structural information about proteins and other biomolecules, and are widely used in Drug Discovery. However, sorting through large numbers of PDB files to find the desired structures can be a time-consuming task. PDB-CAT is a Jupyter Notebook that aims to simplify this process by automatically categorizing the structures based on the type of interaction between atoms in the protein and the ligand, and checking for any mutations in the sequence. PDB-CAT is a program that classifies a group of protein structures into three categories: covalent bonded, non-covalent bonded, and no-bonded. Before classification, the program verifies if there are any mutations in the protein sequence by searching for the amino acid sequence located in the SEQRES lines of the PDB file and comparing it to a reference sequence. The program outputs a file with structures that meet the specified criteria, which can then be used for further analysis, such as virtual screening or molecular docking. The program is easy to use and can be customized to fit the user's specific needs.

## REQUIREMENTS
This notebook uses Python 3.10.9
This program requires the following packages:
- Biopython (1.81)
- os
- re
- csv

You can install these packages using pip:
    $ pip install <library>

*You do not need extensive knowledge of these libraries.*

## USAGE
This program takes PDB files locally saved as an input. It classifies the batch structures based on their type of complexes. The classification is based on the type of interaction between atoms in the structure, it has three categories:

1. Covalent: if the structure contains a LINK information line in the PDB file
2. Non-covalent: if the structure contains at least one non-covalent bond with a peptide
3. No-bond: if the strucutre has not bond at all, or it bonds solvents or ion molecules (described in blacklist)

In this repository, a text file is attached. The text file named *blacklist* and is used to keep track of all groups identified as non-ligand in a protein structure. HET group refers to a non-standard or "heterogeneous" group of atoms within a protein or nucleic acid structure that are not part of the standard building blocks of amino acids or nucleotides. These can include small molecules such as cofactors, and inhibitors that bind to the protein, as well as metal ions and solvent molecules that are present in the crystal structure.

Additionally, before the classification the program checks if there are any mutations in the sequence. To do this, the program searches for the aminoacids sequence (located in SEQRES lines of PDB file) and compares them to a reference sequence. If there are any differences, the program reports that there is a mutation, where is it, how many gaps has the chain sequence and what is the % of identity to reference sequence.

## MODIFICATION

Every code cell where user needs to modify information is marked in 4 stages. In the stage 3, is recommended but not mandatory to check first CSV output. Then, the program can carry on the subsequent functions.

![Image URL](image_documentation/stage_3.png)

## INPUT

The first step before running this notebook is to search the target in https://www.rcsb.org/ webpage. It is needed to download the different structures that match our search and save them in a folder inside the repository directory. Thus to start running the code, the input files are the PDB files in the directory and the *blacklist.txt*. It is necessary to ensure that the blacklist.txt file is located in the working directory before running the program.

## OUTPUT
This code consists of two main parts.

- In the first part, the code reads PDB files and extracts information to create three CSV files:

*sequence_SEQRES*: This file contains the amino acid sequence of the protein from the SEQRES records of the PDB file.

*no_mutated_sequence*: This file contains essential information about the quality of the PDB structure sequence, including the sequence length, gaps with respect to a reference structure, and the % of identity towards a reference structure.

*mutated_sequence*: This file contains essential information about the quality of the PDB structure sequence, including the sequence length, gaps with respect to a reference structure, the percentage of identity towards a reference structure, and residues that do not correspond to the 20 essential amino acids.

- In the second part, the code creates four output folders:

*structures_for_docking*: This folder contains all the PDB files that are classified in *no_mutated_sequence*.

*Covalent_Folder*: This folder contains PDB files of covalently bound protein-ligand complexes.

*Non-Covalent_Folder*: This folder contains PDB files of non-covalently bound protein-ligand complexes.

*No-Bond_Folder*: This folder contains PDB files of unbound proteins (also known as apo-structures).

In addition to the output folders, the code also generates three CSV files:

*Data_Covalent*: This file contains information about covalent link between the protein and the ligand, including the information line (LINK line) in the PDB files.

*Data_No-Bond*: This file contains information about non-covalent bond between the protein and the ligand, including the information line (HET line) in the PDB files.

*Data_Non-Covalent*: This file contains information about non-covalent bond between the protein and ions, solvents or pseudopeptides, including the information line (HET line) in the PDB files.

Overall, this code is designed to analyze and categorize protein-ligand complexes based on their structural and chemical features, providing useful information for drug discovery and other applications in structural biology.

The program outputs the following information:
![Image URL](image_documentation/Scheme.png)

## EXAMPLE
A research group needs to study PLpro *in silico*. In order to do a virtual screening, on or two structures have to be chosen. For this reason, this program will help on a database inspection and classification. 

Previously to run the code, a batch of structures should be downloaded. The users of this program should be familiarized with the target, in this example, PLpro and also, how ligand acts in catalytic site. 

![Image URL](image_documentation/PDB_webpage.png)

*A PDB search is been done, therefore, next step should be *download all*.*

First of all, *in the first stage*(1), the user needs to modify the code, specifically: **path** and **directory**. 

![Image URL](image_documentation/stage_1.png)

Path: where ouput files are created
Directory: where input files should be saved before running the program

In the *second stage*(2), structure and chain of reference should be chosen. In this case, the group choses 7JIW, as it has no mutation, it is in a non-covalent complex and also has a reference length (not short, not large) with no gaps in between.

![Image URL](image_documentation/stage_2.png)

Here the reference structure would be 7JIW.ent from the PDB, and the reference chain would be the A chain.

In the *third stage (3)*, two CSV files have been created. You can go to your path directory and check those sequences in *mutated_sequence* and in *no_mutated_sequence*, or continue executing the code. It should be noted that only PDB files from good_structures will be classified.

Note: "path" refers to the location of the folder on your system where the mentioned CSV files were saved.

For example, the idea is to eliminate rows in *no_mutated_sequence* that have short length or to add some rows of *mutated_sequence* row as it maybe have a mutation isolated from the binding site.

![Image URL](image_documentation/example_modify_first_input.png)

*This is how the first input looks like. From now on, the program will take this information to make the complexes classification. Remember to save as csv!*

Finally, *the forth stage*(4) is to chose a residue of key importance, and also the one that bind with covalent bond. 

![Image URL](image_documentation/stage_4.png)

In this case, it would be the Cys111. 

After the classification:
We see three folders: Covalent_Folder(4), Non-Covalent_Folder(16) and No-Bond_Folder(3); and also, three csv with the binding information: Data_Covalent.csv, Non-Covalent_Folder.csv, Data_No-Bond.csv.

![Image URL](image_documentation/output.png)

## LIMITATIONS

One potential limitation of the PDB-CAT program is that the program requires the user to have a basic understanding of Python programming and the Biopython package, which may be a barrier for some users who are not familiar with these tools.