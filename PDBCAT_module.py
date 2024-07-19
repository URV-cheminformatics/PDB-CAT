
'''
Import libraries
'''

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
    """
    Read a blacklist file and return a list of entries.

    Parameters:
    - blacklist_path (str): Path to the blacklist file.

    Returns:
    - list: List of entries read from the blacklist file.
    """

    with open(blacklist_path, 'r') as f:
        lines = f.readlines()
        blacklist_dict = {}

        for line in lines:
            if not line.startswith('#'):
                columns = line.split()
                key = columns[0].rstrip(',')
                value = columns[1:len(columns)]
                value = " ".join(value)
                blacklist_dict[key.upper()] = value

        blacklist = list(blacklist_dict.keys())  

    return blacklist, blacklist_dict



def search_asym_id(cif_data, entity_id):
    """ 
    From _entity_id we search the _struct_asym.id.
    The _struct_asym.id is also the chain or subunit.
    """

    asym_id = []
    if '_struct_asym' in cif_data:
        for r in cif_data._struct_asym.search('entity_id',entity_id).values():
            asym_id.append(r['id'])
        asym_id = ','.join(asym_id)

    return asym_id



def search_peptide_code(cif_data, asym_id):
    """
    Data items in the PDBX_MOLECULE category identify reference molecules within a PDB entry.
    The value of _pdbx_molecule.prd_id is the PDB accession code for this reference molecule.
    """

    pep_code = ''
    if '_pdbx_molecule' in cif_data:
        for r in cif_data._pdbx_molecule.search('asym_id',asym_id).values():
            pep_code = r['prd_id']
        
        return pep_code
    
    

def search_ligand_code(cif_data, entity_id):
    """
    From the entity_id we search _pdbx_entity_nonpoly.comp_id which is a pointer to the _chem_comp.id
    """

    comp_id = ''
    if '_pdbx_entity_nonpoly' in cif_data:
        for r in cif_data._pdbx_entity_nonpoly.search('entity_id',entity_id).values():
            comp_id = r['comp_id']
        
        return comp_id
    

        
def check_polypeptide(cif_data, entity_id):
    """
    Check if the polymer is a 'polypeptide(L)' or 'polypeptide(D)'. Only polypeptides are considered as the main polymers
    Data items in the ENTITY_POLY category record details about the polymer, such as the type of the polymer, 
    the number of monomers and whether it has nonstandard features.

    _entity_poly.type contains the type of the polymer. It can be: cyclic-pseudo-peptide, other, peptide nucleic acid, 
    polydeoxyribonucleotide, polydeoxyribonucleotide/polyribonucleotide hybrid, polypeptide(D), polypeptide(L), polyribonucleotide,
    polysaccharide(D), polysaccharide(L).
    """

    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id', entity_id).values():
            if r['type'] == 'polypeptide(L)' or r['type'] == 'polypeptide(D)':
                
                return 'Yes'
        
        return 'No'
            


def check_modified_residues(cif_data, asym_id, aminoacid, number):
    """
    Check modified residues in sequence
    """

    modification = ''
    if '_pdbx_struct_mod_residue' in cif_data:
        for r in cif_data._pdbx_struct_mod_residue.search('label_asym_id',asym_id).values():
            if aminoacid.upper() == r['label_comp_id'].upper() and int(number) == int(r['label_seq_id']):
                modification = r['details']
                
                return modification
            


def check_covalent(cif_data, asym_id, chains):
    """
    Check if the ligand has a covalent bond with the protein.
    """

    modified_res_message = ''
    out = ''

    if '_struct_conn' in cif_data:
        for chain in chains.split(','):
            for r in cif_data._struct_conn.search('conn_type_id', 'covale').values():
                if (r.get('ptnr2_label_asym_id') == asym_id and r.get('ptnr1_label_asym_id') == chain) or \
                   (r.get('ptnr1_label_asym_id') == asym_id and r.get('ptnr2_label_asym_id') == chain):
                    partner_key = 'ptnr1' if r.get('ptnr2_label_asym_id') == asym_id else 'ptnr2'
                    comp_id_key = f'{partner_key}_label_comp_id'
                    auth_seq_id_key = f'{partner_key}_auth_seq_id'

                    if 'pdbx_role' in r and r['pdbx_role'][0] == '?':
                        modified_res_message = check_modified_residues(
                            cif_data, chain, r.get(comp_id_key), r.get(auth_seq_id_key))
                    else:
                        modified_res_message = r.get('pdbx_role', 'Unknown role')

                    out = f"{r.get(comp_id_key, 'UnknownComp')}{r.get(auth_seq_id_key, 'UnknownSeq')} {chain}"
                    return out, modified_res_message

    return out, modified_res_message



def search_peptide_seq_one_letter_code(cif_data, entity_id):
    """
    Extract protein sequence
    """

    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            seq_one_letter_code = r['pdbx_seq_one_letter_code'] 
            
        return seq_one_letter_code
    

       
def search_peptide_seq_one_letter_code_can(cif_data, entity_id):
    """
    Search standard-residue one letter sequence
    """

    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            seq_one_letter_code_can = ''.join(r['pdbx_seq_one_letter_code_can'].splitlines())
            
        return seq_one_letter_code_can
    

    
def search_comp_type(cif_data,id):
    """
    For standard polymer components, the type of the monomer.
    Note that monomers that will form polymers are of three types:
    linking monomers, monomers with some type of N-terminal (or 5')
    cap and monomers with some type of C-terminal (or 3') cap.
    """

    if '_chem_comp' in cif_data:
        for r in cif_data._chem_comp.search('id',id).values():
            comp_type = r['type']    
        
        return comp_type   
     


def ligand_info(cif_data,entity_id,asym_id,entity_type,pdbx_description):
    """
    Search peptide ligand info
    """

    ligand_id = ''
    ligand_function = ''
    ligand_type = ''

    pep_code = search_peptide_code(cif_data,asym_id) # Search if the ligand is at The Biologically Interesting Molecule Reference Dictionary (BIRD)
    if pep_code:
        ligand_id = pep_code
        
        for r in cif_data._pdbx_molecule_features.search('prd_id',ligand_id).values():
            ligand_function = r['class']
            ligand_type = r['type']       
   
    elif entity_type == 'polymer':
        seq_one_letter_code = search_peptide_seq_one_letter_code(cif_data,entity_id)
        
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            ligand_type = r['type']
        if seq_one_letter_code:
            ligand_id = seq_one_letter_code
        else:
            ligand_id = pdbx_description    
    
    return ligand_id, ligand_function, ligand_type



def search_branched(cif_data, entity_id):
    """
    Search for branched chain carbohydrates.
    """

    ligand_id = ''
    ligand_type = ''

    if '_entity' in cif_data:
        for r in cif_data._entity.search('id',entity_id).values():    
            if r['type'] == 'branched':
                    ligand_id = (r['pdbx_description']) + ',' + ligand_id
    
    if ligand_id:
        for r in cif_data._pdbx_entity_branch.search('entity_id',entity_id).values():
            ligand_type = (r['type']) + ',' + ligand_type
            
    return ligand_id, ligand_type



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
            n = n+1
            if i != j and j != '-':
                mismatch_location.append(i+str(n)+j)

    return mismatches, mismatch_location, identity, gaps



def process_cif_file(file_path, mutation, blacklist, seq_ref, res_threshold):
    """ 
    Process the cif file analyzing specific data as: ID, chain, number of residues, if it is a complex, 
    if it is peptide-like complex, the type of the bond and the name of the ligand
    """
     
    # Function to process a CIF file and extract the required information
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output = 'cif_wrapper', ignore = ['_atom_site']) # if we are interested in data annotations only and not in coordinates
    cif_data = list(cif_obj.values())[0]

    # Variable initialization
    ligands = []
    ligand_names = []
    ligand_types = []
    ligand_functions = []
    ligand_covalents = []
    ligand_covalents_bond = []
    protein_description = []
    chain = '' 
    num_res = []  
    other_entites = []
    mismatches_list = []
    mismatch_location_list = []
    identity_list  =  [] 
    gaps_list  =  []
    blacklist_id  = []
    ligand_details = {}
    pdb_id = file_path[-8:-4]
    title = cif_data['_struct']['title'][0]

    for i in list(zip(cif_data['_entity']['id'], cif_data['_entity']['type'], cif_data['_entity']['src_method'], cif_data['_entity']['pdbx_description'])):
        entity_id = i[0]
        entity_type = i[1]
        entity_src_method = i[2]
        entity_pdbx_description = i[3]
        asym_id = search_asym_id(cif_data,entity_id)
        
        if entity_type == 'polymer': # if _entity.type is 'polymer' we deduce it is the protein or a peptide ligand  

            if check_polypeptide(cif_data, entity_id) == 'Yes': # 'polypeptide(L)' polymers could contain the protein or peptide ligands
                if entity_src_method == 'syn': # if _entity.src_method is 'syn', we classify synthetic polymers as peptide ligands                                      
                    ligand_id, ligand_function, ligand_type = ligand_info(cif_data, entity_id, asym_id[0], entity_type, entity_pdbx_description)
                    ligands.append(ligand_id)
                    if ligand_id not in ligand_details:
                        ligand_details[ligand_id] = {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                        ligand_details[ligand_id]['type'] = ligand_type
                        ligand_details[ligand_id]['function'] = ligand_function
                        ligand_details[ligand_id]['name'] = entity_pdbx_description
                        ligand_details[ligand_id]['entity'] = entity_id             
                
                else: 
                    seq = search_peptide_seq_one_letter_code_can(cif_data, entity_id)
                    
                    if len(seq) < res_threshold: # Small peptides are considered ligands
                        ligand_id, ligand_function, ligand_type = ligand_info(cif_data, entity_id, asym_id[0], entity_type, entity_pdbx_description)
                        ligands.append(ligand_id)
                        if ligand_id not in ligand_details:
                            ligand_details[ligand_id] = {}
                            ligand_details[ligand_id]['type'] = ligand_type
                            ligand_details[ligand_id]['function'] = ligand_function  
                            ligand_details[ligand_id]['name'] = entity_pdbx_description
                            ligand_details[ligand_id]['entity'] = entity_id                                   
                    else:
                        if chain == '': # chain is the asym_id of the polymers that are not ligands
                            chain = asym_id
                        else:
                            chain = chain + ',' + asym_id
                        
                        protein_description.append(entity_pdbx_description)
                        num_res.append(str(len(seq)))
                        
                        mismatches, mismatch_location, identity, gaps = (None, None, None, None)
                        if mutation == True:
                            mismatches, mismatch_location, identity, gaps = search_for_mutation(seq_ref,seq)

                        mismatches_list.append(mismatches)
                        mismatch_location_list.append(mismatch_location)
                        identity_list.append(identity)
                        gaps_list.append(gaps)  
        
            else:
                j = int(entity_id) - 1
                other_entites.append(cif_data['_entity_poly']['type'][j])

        elif entity_type == 'non-polymer': # non-polymer could be ligands or small-molecules from the medium
            ligand_code = search_ligand_code(cif_data, entity_id) # From the entity_id, search the ligand_code or sequence (for peptides without a code)

            if ligand_code not in blacklist: # Only molecules not present in the blacklist file are considered as ligands
                ligands.append(ligand_code)
                if ligand_code not in ligand_details:
                    ligand_details[ligand_code] = {}
                    comp_type = search_comp_type(cif_data, ligand_code)
                    if comp_type: # If the _chem_comp.type exist, use the _chem_comp.type
                        ligand_details[ligand_code]['type'] = comp_type
                    else:
                        ligand_details[ligand_code]['type'] = 'small-molecule'
                    
                    ligand_details[ligand_code]['function'] = ''  
                    ligand_details[ligand_code]['name'] = entity_pdbx_description
                    ligand_details[ligand_code]['entity'] = entity_id
            
            if ligand_code in blacklist:
                blacklist_id.append(ligand_code)

        if entity_type == 'branched':
            pep_code = search_peptide_code(cif_data,asym_id[0])
            
            if pep_code:
                ligand_id, ligand_function, ligand_type = ligand_info(cif_data,entity_id,asym_id[0],entity_type,entity_pdbx_description)
                ligands.append(ligand_id)
                
                if ligand_id not in ligand_details:
                    ligand_details[ligand_id] = {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                    ligand_details[ligand_id]['type'] = ligand_type
                    ligand_details[ligand_id]['function'] = ligand_function
                    ligand_details[ligand_id]['name'] = entity_pdbx_description
                    ligand_details[ligand_id]['entity'] = entity_id      

            else:
                ligand_id, ligand_type = search_branched(cif_data,entity_id) 
                ligands.append(ligand_id)
                ligand_function = ''
                
                if ligand_id not in ligand_details:
                    ligand_details[ligand_id] = {}       
                    ligand_details[ligand_id]['type'] = ligand_type
                    ligand_details[ligand_id]['name'] = ligand_id
                    ligand_details[ligand_id]['function'] = ligand_function
                    ligand_details[ligand_id]['entity'] = entity_id

    nchain = chain.count(',') + 1
    num_res = ', '.join(num_res)
    protein_description = '\n'.join(protein_description)
    blacklist_id = '\n'.join(blacklist_id)
    other_entites = '\n'.join(other_entites)

    # Search if ligands are covalent bonded
    for i in ligands: 
        asym_id = search_asym_id(cif_data,ligand_details[i]['entity'])
        covalent_bond = ''
        covalent = ''
        bond = ''
        message = ''

        for id_a in asym_id.split(','): # If an entity has several chains, only one chain needs to have a covalent bond to be considered a covalent bond
            covalent_bond,message = check_covalent(cif_data,id_a,chain)
            if covalent_bond:
                if message:
                    covalent = '(' + str(message) + ')'
                    if message == "Unknown role":
                        covalent = 'Yes'
                else:
                    covalent = 'Yes'
                
                bond = bond + '; ' + covalent_bond
            
        if covalent:
            ligand_details[i]['covalent'] = covalent
            ligand_details[i]['covalent_bond'] = bond[2:] 

        else:
            ligand_details[i]['covalent'] = 'No'
            ligand_details[i]['covalent_bond'] = ''
               

    for i in ligand_details.keys():
        ligand_names.append(ligand_details[i]['name'])
        ligand_types.append(ligand_details[i]['type'])
        ligand_functions.append(ligand_details[i]['function'])
        ligand_covalents.append(ligand_details[i]['covalent'])
        ligand_covalents_bond.append(ligand_details[i]['covalent_bond'])

    ligand_names = '\n'.join(ligand_names)
    ligand_types = '\n'.join(ligand_types)
    ligand_functions = '\n'.join(ligand_functions)
    ligand_covalents  =  '\n'.join(ligand_covalents)
    ligand_covalents_bond  =  '\n'.join(ligand_covalents_bond)
    ligands = '\n'.join(ligands)
  
    if not ligands: complex = 'No'
    else: complex = 'Yes'

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
        "Title": title,
        "Protein_description": protein_description,
        "NChain": nchain,
        "Chain_ID": chain,
        "Num_Res": num_res,
        "Complex": complex,
        "Discarted_Ligands": blacklist_id,
        "Other_Entities": other_entites,
        "Ligand": ligands,
        "Ligand_names": ligand_names,
        "Ligand_types": ligand_types,
        "Ligand_functions": ligand_functions,
        "Covalent_Bond": ligand_covalents,
        "Bond": ligand_covalents_bond,
        "Mutation": mismatches,
        "Mutation_Location": mismatch_location,
        "Identity": identity,
        "Gaps": gaps
    }



def mutation_classification(directory_path, information_df, output_path):
    """
    Classification based on the information contained in a dataframe. 
    Downloading PDB files into different folders depend whether there is a mutation.
    """
       
    mutated_list = []
    no_mutated_list = []

    current_datetime = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

    with open(information_df, 'r') as f:
        df = pd.read_csv(f, dtype = {'Mutation':str})
        for index, row in df.iterrows():
            if row['Mutation'] == '0':
                no_mutated_list.append(row['PDB_ID'])
            else:
                mutated_list.append(row['PDB_ID'])

    origin_path = directory_path

    mut_path = output_path + '/Mutated/'
    if os.path.exists(mut_path):
        mut_path = output_path + '/Mutated_' + current_datetime + '/'
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
        pdb_id = mutated
        shutil.copy(origin_path + (pdb_id+".cif"), mut_path+(pdb_id + ".cif"))
        
    for no_mutated in no_mutated_list:
        pdb_id = no_mutated
        shutil.copy(origin_path + (pdb_id+".cif"), non_mut_path + (pdb_id +".cif"))
    
    return no_mutated_list, non_mut_path 



def bond_classification(directory_path, information_df, no_mutated_list, output_path, mutation):
    """
    This function classifies between free enzymes and complexes.
    And for each complex, it classifies depending on their bonds: covalent or non-covalent.
    """

    free_enzyme_list = []
    covalent_list = []
    non_covalent_list = []
    type_bond = ''

    with open(information_df, 'r') as f:
        df = pd.read_csv(f)
        for index, row in df.iterrows():
            if row['Complex'] == 'No':
                if row['PDB_ID'] in no_mutated_list:
                    free_enzyme_list.append(row['PDB_ID'])
                    
            else:
                if row['PDB_ID'] in no_mutated_list:
                    type_bond = row['Covalent_Bond']
                    if 'Yes' in type_bond:
                        covalent_list.append(row['PDB_ID'])
                    else: 
                        non_covalent_list.append(row['PDB_ID'])

    if mutation == False:
        origin_path = directory_path
    elif mutation == True:
        origin_path = output_path

    free_path = output_path  + 'APO/'
    covalent_path = output_path  + 'Covalent/'
    non_covalent_path = output_path + 'Non-Covalent/'

    # Create directories if they don't exist
    for path in [free_path, covalent_path, non_covalent_path]:
        if not os.path.exists(path):
            os.makedirs(path)

    for free in free_enzyme_list:
        pdb_id = free + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), free_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), free_path+(pdb_id))

    for covalent in covalent_list:
        pdb_id = covalent + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), covalent_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), covalent_path+(pdb_id))

    for non_covalent in non_covalent_list:
        pdb_id = non_covalent + ".cif"
        if mutation == True:
            shutil.move(origin_path+(pdb_id), non_covalent_path+(pdb_id))
        if mutation == False:
            shutil.copy(origin_path+(pdb_id), non_covalent_path+(pdb_id))
    
    return
