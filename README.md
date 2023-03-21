# PROTIC

PROTIC is a Python program that identifies no-mutated sequences of a specific target by downloading and classifying PDB files. The program categorizes the batch structures based on the type of interaction between atoms in the protein and the ligand into three categories: Covalent, Non-covalent, and Free enzymes. The program checks for any mutations in the sequence before classification by searching for the amino acid sequence (located in SEQRES lines of PDB file) and comparing them to a reference sequence.

## REQUIREMENTS
This program requires the following packages:
- Biopython
- os
- re
- csv

You can install these packages using pip:
    $ pip install <library>

## USAGE
This program takes PDB files locally saved as an input. It classifies the batch structures based on their type of complexes. The classification is based on the type of interaction between atoms in the structure, it has three categories:

1. Covalent: if the structure contains only covalent bonds
2. Non-covalent: if the structure contains only non-covalent bonds
3. Free: if the strucutre does not bond to any ligand

Additionally, before the classification the program checks if there are any mutations in the sequence. To do this, the program searches for the aminoacids sequence (located in SEQRES lines of PDB file) and compares them to a reference sequence. If there are any differences, the program reports that there is a mutation, where is it, how many gaps has the chain sequence and what is the % of identity to reference sequence.

## MODIFICATION

Every code cell where user needs to modify information is marked with numeric point. There are 4 points. In point 3, user need to make a STOP and to modify the first output. Then, the program can carry on the previous functions.

## INPUT
This code has two parts. 

In the first part, three CSV should be created. 
1. sequence_SEQRES
2. good_structures_info
3. mutation_info

Information about mutations, gaps and %identity.

In the second part, the programs create 4 folders:
1. structures_for_docking
2. Covalent_Folder
3. NonCovalent_Folder
4. Free_Folder

PDB files classify for complex type.

and also 3 CSV:
1. DataCovalent
2. DataFree
3. DataNonCovalent

Information about bond, where is the bond and to how it is made.

## OUTPUT
The program outputs the following information:

The classification of the structure (covalent, non-covalent, or no bond).
If the structure is classified as non-covlent, the program also reports the type of non-covalent interaction present in the structure (hydrogen bond, salt bridge, etc). 
If the program detects a mutatiion, it reports the position of the mutation and the aminoacids change, gaps and identity %. Also, you can access to a csv with all no-mutated structures to check the length and gaps these.

## EXAMPLE
A research group needs to study PLpro *in silico*. In order to do a virtual screening, on or two structures have to be chosen. For this reason, this program will help on a database inspection and classification. 

Previously to run the code, a batch of structures should be downloaded. The users of this program should be familiarized with the target, in this example, PLpro and also, how ligand acts in catalytic site. 

![Image URL](/home/ariadna/PROTIC/image_documentation/PDB_webpage.png)
A PDB search is been done, therefore, next step should be *download all*.

First of all, in the first point, the user needs to modify the code, specifically: **path** and **directory**. 
Path: where ouput files are created
Directory: where input files should be saved before running the program

In the second point, structure and chain of reference should be chosen. In this case, the group choses 7JIW, as it has no mutation, it is in a non-covalent complex and also has a reference length (not short, not large) with no gaps in between.

![Image URL](/home/ariadna/PROTIC/image_documentation/7jiw_assembly-1.jpeg)
7JIW 3D structure in PDB

In third point, we are going to **STOP** as two csv files has been created. You should go to your **path** directory and check those sequences in *mutation_info* and in *good_structures*. It has to be known that only *good_structures* PDB files will be classify.

For example, the idea is to eliminate *good_structures* with short length or to add some *mutation_info* row as it maybe have a mutation isolated from the binding site.

Finally, the forth point is to chose a residue of key importance, and also the one that bind with covalent bond. In this case, it is Cys111.

After the classification:
We see three folders: covalent(2), non_covalent(15) and free(3); and also, three csv with the binding information: DataCovalent, DataNonCovalent, DataFree.

