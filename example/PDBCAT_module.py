# PDBCATmodule

'''
Import libraries
'''
from PDBCAT_module import *
from pdbecif.mmcif_io import CifFileReader
from pdbecif.mmcif_tools import MMCIF2Dict
import pandas as pd
import time
import re
import os
import shutil
from Bio.Align import PairwiseAligner 
from Bio.PDB import *  
from datetime import datetime

'''
PDB-CAT functions
'''

def read_blacklist(blacklist_path):

    with open(blacklist_path, 'r') as f:
        linea = f.read()
    linea=linea.replace(' ','')
    blacklist = linea.strip('\n').split(',')

    return blacklist


def search_asym_id(cif_data,entity_id):
    """ 
    From _entity_id we search the _struct_asym.id.
    The _struct_asym.id is also the chain or subunit.
    """

    asym_id=[]
    if '_struct_asym' in cif_data:
        for r in cif_data._struct_asym.search('entity_id',entity_id).values():
            asym_id.append(r['id'])
        asym_id=','.join(asym_id)

    return asym_id


def search_peptide_code(cif_data,asym_id):
    """
    Data items in the PDBX_MOLECULE category identify reference molecules within a PDB entry.
    The value of _pdbx_molecule.prd_id is the PDB accession code for this reference molecule.
    """

    pep_code=''
    if '_pdbx_molecule' in cif_data:
        for r in cif_data._pdbx_molecule.search('asym_id',asym_id).values():
            pep_code=r['prd_id']
        return pep_code
    

def search_ligand_code(cif_data, entity_id):
    """
    From the entity_id we search _pdbx_entity_nonpoly.comp_id which is a pointer to the _chem_comp.id
    """

    comp_id=''
    if '_pdbx_entity_nonpoly' in cif_data:
        for r in cif_data._pdbx_entity_nonpoly.search('entity_id',entity_id).values():
            comp_id=r['comp_id']
        return comp_id
    
    
def check_polypeptide(cif_data, entity_id):
    """
    Check if the polymer is a 'polypeptide(L)'. Only polypeptides are considered as the main polymers
    Data items in the ENTITY_POLY category record details about the polymer, such as the type of the polymer, 
    the number of monomers and whether it has nonstandard features.
    _entity_poly.type contains the type of the polymer. It can be: cyclic-pseudo-peptide, other, peptide nucleic acid, 
    polydeoxyribonucleotide, polydeoxyribonucleotide/polyribonucleotide hybrid, polypeptide(D), polypeptide(L), polyribonucleotide,
    polysaccharide(D), polysaccharide(L).
    """

    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            if r['type'] == 'polypeptide(L)':
                return 'Yes'
            else:
                return 'No'
            
            
def check_covalent(cif_data,asym_id): 
    """
    Check if the ligand has a covalent bond with the protein
    """

    if '_struct_conn' in cif_data:
        for r in cif_data._struct_conn.search('conn_type_id','covale').values(): 
            if r['ptnr2_label_asym_id'] == asym_id and r['ptnr2_label_asym_id'] != r['ptnr1_label_asym_id']: 
                return r['ptnr1_label_comp_id']+r['ptnr1_auth_seq_id'] 
            if r['ptnr1_label_asym_id'] == asym_id and r['ptnr2_label_asym_id'] != r['ptnr1_label_asym_id']:
                return r['ptnr2_label_comp_id']+r['ptnr2_auth_seq_id'] 

def search_peptide_seq_one_letter_code(cif_data,entity_id):
    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            seq_one_letter_code = r['pdbx_seq_one_letter_code'] 
            
        return seq_one_letter_code
    
    
def search_peptide_seq_one_letter_code_can(cif_data,entity_id):
    """
    Search standard-residue one letter sequence
    """
    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            seq_one_letter_code_can = ''.join(r['pdbx_seq_one_letter_code_can'].splitlines())
            
        return seq_one_letter_code_can
    

def search_for_mutation(reference_seq, one_letter_seq):
    """ 
    Search for mutations, where is the mismatch, how many gaps there are and the identity of the sequence
    """

    # Define aligner variables from PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 5
    aligner.mismatch_score = -3
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -2
    
    # Perform pairwise alignment
    alignments = aligner.align(reference_seq, one_letter_seq)
    alignment = alignments[0] # We select the first alignment. It could be more than 1 aligment with same score 
    identity = '{:.2f}'.format(alignment.counts()[1]*100/(len(alignment[0, :])))
    gaps = alignment.counts()[0]

    # Extract the count of mismatches
    mismatches = alignment.counts()[2]
    # Where is the mismatch
    n = 0
    mismatch_location = []

    for i, j in zip(alignment[0, :], alignment[1, :]):
        if i != '-':
            n=n+1
            if i!=j and j!='-':
                mismatch_location.append(i+str(n)+j)
    
    return mismatches, mismatch_location, identity, gaps


def process_cif_file(file_path, mutation, blacklist, seq_ref, res_threshold):
    """ 
    Process the cif file analyzing specific data as: ID, chain, number of residues, if it is a complex, 
    if it is peptide-like complex, the type of the bond and the name of the ligand
    """

    # Process CIF file
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output='cif_wrapper', ignore=['_atom_site']) # if we are interested in data annotations only and not in coordinates, we can discard the category _atom_site with ignore=['_atom_site'].
    cif_data = list(cif_obj.values())[0]

    # Variable initialization
    ligands=[]
    ligand_names=[]
    pols=[]
    chain=[] # Only of the polymer proteins
    num_res=[]  
    mismatches_list = []
    mismatch_location_list = []
    identity_list = [] 
    gaps_list = []
    blacklist_id =[]
    peptide_like='No' # if a complex contains one peptide_like ligand, this variable will be 'Yes'
    covalent='No' # if one ligand is found to have a covalent bond with the protein, the pdb id is considered covalent='Yes'
    covalent_bond=''
    bond=''

    pdb_id = cif_data['_entry']['id'][0] # The PDB ID
    title = cif_data['_struct']['title'][0]

    for i in list(zip(cif_data['_entity']['id'],cif_data['_entity']['type'],cif_data['_entity']['src_method'],cif_data['_entity']['pdbx_description'])):
        entity_id = i[0]
        if i[1] == 'polymer': # if _entity.type is 'polymer' we deduce it is the protein or a peptide ligand  

            if i[2] == 'syn': # if _entity.src_method is 'syn', we classify synthetic polymers as peptide ligands
                peptide_like='Yes'
                ligand_names.append(i[3])
                asym_id=search_asym_id(cif_data,entity_id)
                for id_a in asym_id.split(','):
                    covalent_bond=check_covalent(cif_data,id_a)
                    if covalent_bond:
                        covalent='Yes'
                        bond=covalent_bond                    

                    pep_code=search_peptide_code(cif_data,id_a) # Search if the peptide ligand has a prd.id code, such as 'PRD_000338'
                    if pep_code != None:
                        if pep_code not in ligands: ligands.append(pep_code)
                    else: # Ih the peptide code has not a prd.id code, search the _entity_poly.pdbx_seq_one_letter_code
                        seq_one_letter_code=search_peptide_seq_one_letter_code(cif_data,entity_id)
                        if seq_one_letter_code != None:
                            if seq_one_letter_code not in ligands: ligands.append(seq_one_letter_code)
                        else:
                            if i[3] not in ligands: ligands.append(i[3])
                
            else:
                if check_polypeptide(cif_data, entity_id) == 'Yes': # Only 'polypeptide(L)' polymers are considered as the main polymers of the structure
                    seq = search_peptide_seq_one_letter_code_can(cif_data,entity_id)
                      
                    # check if its a peptide 
                    if len(seq) < res_threshold:
                        peptide_like='Yes' 
                        ligand_names.append(i[3])
                        asym_id=search_asym_id(cif_data,entity_id)
                        for id_a in asym_id.split(','):
                            covalent_bond=check_covalent(cif_data,id_a)
                            if covalent_bond:
                                covalent='Yes'
                                bond=covalent_bond 
                            
                            seq_one_letter_code=search_peptide_seq_one_letter_code(cif_data,entity_id)
                            if seq_one_letter_code != None:
                                if seq_one_letter_code not in ligands: ligands.append(seq_one_letter_code)
                            else:
                                if i[3] not in ligands: ligands.append(i[3])
                                
                    else:
                        chain.append(search_asym_id(cif_data,entity_id))  # chain is the asym_id of the polymers that are not ligands
                        pols.append(i[3])
                        num_res.append(str(len(seq)))
                        
                        mismatches, mismatch_location, identity, gaps = (None, None, None, None)
                        if mutation == True:
                            mismatches, mismatch_location, identity, gaps = search_for_mutation(seq_ref,seq) 
                        
                        mismatches_list.append(mismatches)
                        mismatch_location_list.append(mismatch_location)
                        identity_list.append(identity)
                        gaps_list.append(gaps)  

        elif i[1] == 'non-polymer': # non-polymer could be ligands or small-molecules from the medium
            ligand_code=search_ligand_code(cif_data,entity_id) # From the entity_id, search the ligand_code or sequence (for peptides without a code)
            if ligand_code not in blacklist: # Only molecules not present in the blacklist file are considered as ligands
                ligands.append(ligand_code)
                ligand_names.append(i[3])
                asym_id=search_asym_id(cif_data,entity_id)
                for id_a in asym_id.split(','):
                    covalent_bond=check_covalent(cif_data,id_a)
                    if covalent_bond:
                        covalent='Yes'
                        bond=covalent_bond
            if ligand_code in blacklist:
                blacklist_id.append(ligand_code)

    chain=','.join(chain) # We join all the chains
    nchain=chain.count(',')+1
    num_res=', '. join(num_res)
    ligand_names='\n'. join(ligand_names)
    pols='\n'.join(pols)
    blacklist_id='\n'. join(blacklist_id)
    
    if not ligands: complex='No'
    else: complex ='Yes'

    ligands='\n'.join(ligands)
    if mutation == True:
        mismatch_location = ', '.join(mismatch_location)
    else:
        mismatch_location = None
    
    mismatches = '\n'.join(map(str, mismatches_list))
    mismatch_location = '\n'.join(map(str, mismatch_location_list))
    identity = '\n'.join(map(str,identity_list))
    gaps = '\n'.join(map(str,gaps_list))

    return {
        "PDB_ID": pdb_id,
        "Title":title,
        "Protein_description": pols,
        "NChain":nchain,
        "Chain_ID": chain,
        "Num_Res": num_res,
        "Complex": complex,
        "Descarted_Ligands": blacklist_id,
        "Ligand": ligands,
        "Ligand_names": ligand_names,
        "Peptide_like": peptide_like,
        "Covalent_Bond": covalent,
        "Bond": bond,
        "Mutation": mismatches,
        "Mutation_Location": mismatch_location,
        "Identity": identity,
        "Gaps": gaps
    }


def mutation_classification(information_df, output_path):
    """
    Classification based on the information contained in a dataframe. 
    Downloading PDB files into different folders depend whether there is a mutation.
    """
       
    mutated_list = []
    no_mutated_list = []

    current_datetime = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

    with open(information_df, 'r') as f:
        df = pd.read_csv(f)
        for index, row in df.iterrows():
            if row['Mutation'] == '0':
                no_mutated_list.append(row['PDB_ID'])
            else:
                mutated_list.append(row['PDB_ID'])

    origin_path = os.getcwd() +'/cif/'

    mut_path = output_path + '/Mutated/'
    if os.path.exists(mut_path):
        mut_path = output_path + '/Mutated_' + current_datetime +'/'
        os.makedirs(mut_path)
    else: 
        os.makedirs(mut_path)

    non_mut_path = output_path + '/Non-Mutated/'
    if os.path.exists(non_mut_path):
        non_mut_path = output_path + '/Non-Mutated_' + current_datetime + '/'
        os.makedirs(non_mut_path)
    else:
        os.makedirs(non_mut_path)


    for mutated in mutated_list:
        pdb_id = mutated.lower()
        shutil.copy(origin_path+(pdb_id+".cif"), mut_path+(pdb_id+".cif"))

    for no_mutated in no_mutated_list:
        pdb_id = no_mutated.lower()
        shutil.copy(origin_path+(pdb_id+".cif"), non_mut_path +(pdb_id+".cif"))


    return no_mutated_list, non_mut_path 


def bond_classification(information_df, no_mutated_list, output_path, mutation):
    """
    This function classifies between free enzymes and complexes.
    And for each complex, it classifies depending on their bonds: covalent or non-covalent.
    """

    free_enzyme_list = []
    covalent_list = []
    non_covalent_list = []

    with open(information_df, 'r') as f:
        df = pd.read_csv(f)
        for index, row in df.iterrows():
            if row['Complex'] == 'No':
                if row['PDB_ID'] in no_mutated_list:
                    free_enzyme_list.append(row['PDB_ID'])
                    
            else:
                if row['PDB_ID'] in no_mutated_list:
                    if row['Covalent_Bond'] == 'Yes':
                        covalent_list.append(row['PDB_ID'])
                    elif row['Covalent_Bond'] == 'No':
                        non_covalent_list.append(row['PDB_ID'])

    if mutation == False:
        origin_path = os.getcwd() +'/cif/'
    elif mutation == True:
        origin_path = output_path

    free_path = output_path  + 'Free/'
    covalent_path = output_path  + 'Covalent/'
    non_covalent_path = output_path + 'Non-Covalent/'

    # Create directories if they don't exist
    for path in [free_path, covalent_path, non_covalent_path]:
        if not os.path.exists(path):
            os.makedirs(path)

    for free in free_enzyme_list:
        pdb_id = free.lower() + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), free_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), free_path+(pdb_id))

    for covalent in covalent_list:
        pdb_id = covalent.lower() + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), covalent_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), covalent_path+(pdb_id))

    for non_covalent in non_covalent_list:
        pdb_id = non_covalent.lower() + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), non_covalent_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), non_covalent_path+(pdb_id))


    return