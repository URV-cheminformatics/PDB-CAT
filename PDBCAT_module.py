# PDBCATmodule

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
                value = columns[1]
                blacklist_dict[key.upper()] = value

        blacklist = list(blacklist_dict.keys())  

    return blacklist, blacklist_dict


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
            
def check_modified_residues(cif_data,asym_id,aminoacid,number):
    modification=''
    if '_pdbx_struct_mod_residue' in cif_data:
        for r in cif_data._pdbx_struct_mod_residue.search('label_asym_id',asym_id).values():
            if aminoacid.upper() == r['label_comp_id'].upper() and int(number) == int(r['label_seq_id']):
                modification = r['details']
                return modification

def check_covalent(cif_data,asym_id,chains): 
    """
    Check if the ligand has a covalent bond with the protein
    """
    modified_res_message=''
    if '_struct_conn' in cif_data:
        for chain in chains.split(','):
            for r in cif_data._struct_conn.search('conn_type_id','covale').values():
                if r['ptnr2_label_asym_id'] == asym_id and r['ptnr1_label_asym_id'] == chain: 
                    modified_res_message=check_modified_residues(cif_data,chain,r['ptnr1_label_comp_id'],r['ptnr1_auth_seq_id'])
                    out=r['ptnr1_label_comp_id']+r['ptnr1_auth_seq_id']+ ' '+chain
                    return out,modified_res_message 
                if r['ptnr1_label_asym_id'] == asym_id and r['ptnr2_label_asym_id'] == chain:
                    modified_res_message=check_modified_residues(cif_data,chain,r['ptnr2_label_comp_id'],r['ptnr2_auth_seq_id'])
                    out=r['ptnr2_label_comp_id']+r['ptnr2_auth_seq_id']+' '+chain
                    return out,modified_res_message 
                
    return '',''

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
    
def search_comp_type(cif_data,id):
     if '_chem_comp' in cif_data:
        for r in cif_data._chem_comp.search('id',id).values():
            comp_type = r['type']    
        return comp_type   

def ligand_info(cif_data,entity_id,asym_id,entity_type,pdbx_description):
    """
    Search peptide ligand info
    """

    ligand_id=''
    ligand_function=''
    ligand_type=''

    pep_code=search_peptide_code(cif_data,asym_id) # Search if the ligand is at The Biologically Interesting Molecule Reference Dictionary (BIRD). Then, it has a prd.id code, such as 'PRD_000338'
    if pep_code:
        ligand_id = pep_code
        for r in cif_data._pdbx_molecule_features.search('prd_id',ligand_id).values():
            ligand_function=r['class']
            ligand_type=r['type']       
   
    elif entity_type == 'polymer':
        seq_one_letter_code=search_peptide_seq_one_letter_code(cif_data,entity_id)
        for r in cif_data._entity_poly.search('entity_id',entity_id).values():
            ligand_type = r['type']
        if seq_one_letter_code:
            ligand_id=seq_one_letter_code
        else:
            ligand_id=pdbx_description    
    
    return ligand_id, ligand_function, ligand_type

def search_branched(cif_data, entity_id):
    
    branched_code = ''
    branched_type = ''
    if '_entity' in cif_data:
        for r in cif_data._entity.search('id',entity_id).values():    
            if r['type'] == 'branched':
                    branched_code=r['pdbx_description']
    
    if branched_code:
        for r in cif_data._pdbx_entity_branch.search('entity_id',entity_id).values():
            branched_type = r['type']
            
    return branched_code, branched_type

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

    # Función para procesar un archivo CIF y extraer la información requerida
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output='cif_wrapper', ignore=['_atom_site']) # if we are interested in data annotations only and not in coordinates, we can discard the category _atom_site with ignore=['_atom_site'].
    cif_data = list(cif_obj.values())[0]

    # Variable initialization
    ligands=[]
    ligand_names=[]
    ligand_types=[]
    ligand_functions=[]
    ligand_covalents=[]
    ligand_covalents_bond=[]
    protein_description=[]
    chain='' # Only of the polymer proteins
    num_res=[]  
    mismatches_list = []
    mismatch_location_list = []
    identity_list = [] 
    gaps_list = []
    blacklist_id =[]
    ligand_details={}
    branched_details = {}
    branched = []
    branched_names=[]
    branched_types=[]
    branched_functions=[]
    branched_covalents=[]
    branched_covalents_bond=[]
    pdb_id = cif_data['_entry']['id'][0] # The PDB ID
    title = cif_data['_struct']['title'][0]

    for i in list(zip(cif_data['_entity']['id'],cif_data['_entity']['type'],cif_data['_entity']['src_method'],cif_data['_entity']['pdbx_description'])):
        entity_id = i[0]
        entity_type = i[1]
        entity_src_method = i[2]
        entity_pdbx_description = i[3]
        asym_id=search_asym_id(cif_data,entity_id)
        
        if entity_type == 'polymer': # if _entity.type is 'polymer' we deduce it is the protein or a peptide ligand  
            
            if check_polypeptide(cif_data, entity_id) == 'Yes': # 'polypeptide(L)' polymers could contain the protein or peptide ligands
                if entity_src_method == 'syn': # if _entity.src_method is 'syn', we classify synthetic polymers as peptide ligands                                      
                    ligand_id, ligand_function, ligand_type = ligand_info(cif_data,entity_id,asym_id[0],entity_type,entity_pdbx_description)
                    ligands.append(ligand_id)
                    if ligand_id not in ligand_details:
                        ligand_details[ligand_id]={}
                        ligand_details[ligand_id]['type'] = ligand_type
                        ligand_details[ligand_id]['function'] = ligand_function
                        ligand_details[ligand_id]['name'] = entity_pdbx_description
                        ligand_details[ligand_id]['entity'] = entity_id             
                else: # Only 'polypeptide(L)' polymers, not synthetics are considered as the main polymers of the structure
                    seq = search_peptide_seq_one_letter_code_can(cif_data,entity_id)
                    if len(seq) < res_threshold: # Small peptides are considered ligands
                        ligand_id, ligand_function, ligand_type = ligand_info(cif_data,entity_id,asym_id[0],entity_type,entity_pdbx_description)
                        ligands.append(ligand_id)
                        if ligand_id not in ligand_details:
                            ligand_details[ligand_id]={}
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

        elif entity_type == 'non-polymer': # non-polymer could be ligands or small-molecules from the medium
            
            ligand_code=search_ligand_code(cif_data,entity_id) # From the entity_id, search the ligand_code or sequence (for peptides without a code)
            if ligand_code not in blacklist: # Only molecules not present in the blacklist file are considered as ligands
                ligands.append(ligand_code)
                if ligand_code not in ligand_details:
                    ligand_details[ligand_code]={}
                    comp_type = search_comp_type(cif_data,ligand_code)
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
            
            pep_code = search_peptide_code(cif_data,asym_id)
            if pep_code:
                branched_code = pep_code
                branched.append(branched_code)
                for r in cif_data._pdbx_molecule_features.search('prd_id',branched_code).values():
                    branched_function=r['class']
                    branched_type=r['type'] 

            else:
                branched_code, branched_type = search_branched(cif_data,entity_id) 
                branched.append(branched_code)
                branched_type = ''
                branched_function= ''

            if branched_code not in branched_details:
                branched_details[branched_code]={}           
                branched_details[branched_code]['type'] =  branched_type
                branched_details[branched_code]['name'] = entity_pdbx_description
                branched_details[branched_code]['function'] = branched_function
                branched_details[branched_code]['entity'] = entity_id

    nchain=chain.count(',')+1
    num_res=', '.join(num_res)
    protein_description='\n'.join(protein_description)
    blacklist_id='\n'.join(blacklist_id)
    
    if not ligands: complex='No'
    else: complex ='Yes'

    for i in ligands: # Search if ligands are covalent bonded
        asym_id=search_asym_id(cif_data,ligand_details[i]['entity'])
        covalent_bond=''
        covalent=''
        bond=''
        message=''
        for id_a in asym_id.split(','): # If an entity has several chains, only one chain needs to have a covalent bond to be considered a covalent bond
            covalent_bond,message =check_covalent(cif_data,id_a,chain)
            if covalent_bond:
                if message:
                    covalent = 'Yes'+ ' ('+str(message)+')' # The message can contain for example that it is a glycosylation
                else:
                    covalent='Yes'
                bond=covalent_bond
        if covalent[0:3] == 'Yes':
            ligand_details[i]['covalent']= covalent
            ligand_details[i]['covalent_bond']=bond
        else:
            ligand_details[i]['covalent']='No'
            ligand_details[i]['covalent_bond']=''
             
    ligands='\n'.join(ligands)
    if mutation == True:
        mismatch_location = ', '.join(mismatch_location)
    else:
        mismatch_location = None

    for i in ligand_details.keys():
        ligand_names.append(ligand_details[i]['name'])
        ligand_types.append(ligand_details[i]['type'])
        ligand_functions.append(ligand_details[i]['function'])
        ligand_covalents.append(ligand_details[i]['covalent'])
        ligand_covalents_bond.append(ligand_details[i]['covalent_bond'])

    ligand_names='\n'.join(ligand_names)
    ligand_types='\n'.join(ligand_types)
    ligand_functions='\n'.join(ligand_functions)
    ligand_covalents = '\n'.join(ligand_covalents)
    ligand_covalents_bond = '\n'.join(ligand_covalents_bond)

    for i in branched:
        asym_id=search_asym_id(cif_data,branched_details[i]['entity'])

        covalent_branched=''
        covalent_branched_bond=''
        branched_covalent=''
        message=''
        
        for id_a in asym_id.split(','): 
            covalent_branched, message = check_covalent(cif_data,asym_id,chain)
            if covalent_branched:
                if message:
                    branched_covalent = 'Yes'+ ' ('+str(message)+')'
                else:
                    branched_covalent = "Yes"

                covalent_branched_bond = covalent_branched

            if branched_covalent[0:3] == 'Yes':
                branched_details[i]['covalent']= branched_covalent
                branched_details[i]['covalent_bond']=covalent_branched_bond
            else:
                branched_details[i]['covalent']='No'
                branched_details[i]['covalent_bond']=''

    for i in branched_details.keys():
        branched_names.append(branched_details[i]['name'])
        branched_types.append(branched_details[i]['type'])
        branched_functions.append(branched_details[i]['function'])
        branched_covalents.append(branched_details[i]['covalent'])
        branched_covalents_bond.append(branched_details[i]['covalent_bond'])

    branched = '\n'.join(branched)
    branched_names = '\n'.join(branched_names)
    branched_types = '\n'.join(branched_types)
    branched_functions = '\n'.join(branched_functions)
    branched_covalents = '\n'.join(branched_covalents)
    branched_covalents_bond = '\n'.join(branched_covalents_bond)

    mismatches = '\n'.join(map(str, mismatches_list))
    mismatch_location = '\n'.join(map(str, mismatch_location_list))
    identity = '\n'.join(map(str,identity_list))
    gaps = '\n'.join(map(str,gaps_list))

    return {
        "PDB_ID": pdb_id,
        "Title":title,
        "Protein_description": protein_description,
        "NChain":nchain,
        "Chain_ID": chain,
        "Num_Res": num_res,
        "Complex": complex,
        "Descarted_Ligands": blacklist_id,
        "Branched": branched,
        "Branched_name": branched_names,
        "Branched_type": branched_types,
        "Branched_functions": branched_functions,
        "Branched_Covalent": branched_covalents,
        "Branched_Bond": branched_covalents_bond,
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
        df = pd.read_csv(f)
        for index, row in df.iterrows():
            if row['Mutation'] == '0':
                no_mutated_list.append(row['PDB_ID'])
            else:
                mutated_list.append(row['PDB_ID'])

    origin_path = directory_path

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


def bond_classification(directory_path, information_df, no_mutated_list, output_path, mutation):
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
        origin_path = directory_path
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