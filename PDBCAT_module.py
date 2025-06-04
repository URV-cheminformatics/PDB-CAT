# PDB-CAT program v3.0.0 

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
from Bio import SeqIO
from datetime import datetime
import numpy as np
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import psutil
import threading

# Shared counter for processed files
valid_count = multiprocessing.Value('i', 0)
processed_count = multiprocessing.Value('i', 0)
impression = False

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
    Check if the polymer is a 'polypeptide(L)'. Only polypeptides are considered as the main polymers
    Data items in the ENTITY_POLY category record details about the polymer, such as the type of the polymer, 
    the number of monomers and whether it has nonstandard features.

    _entity_poly.type contains the type of the polymer. It can be: cyclic-pseudo-peptide, other, peptide nucleic acid, 
    polydeoxyribonucleotide, polydeoxyribonucleotide/polyribonucleotide hybrid, polypeptide(D), polypeptide(L), polyribonucleotide,
    polysaccharide(D), polysaccharide(L).
    """

    if '_entity_poly' in cif_data:
        for r in cif_data._entity_poly.search('entity_id', entity_id).values():
            if r['type'] == 'polypeptide(L)':
                return 'Yes'
            else:
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
    covalent = ''

    if '_struct_conn' in cif_data:
        for chain in chains.split(','):
            for r in cif_data._struct_conn.search('conn_type_id', 'covale').values():
                if (r.get('ptnr2_label_asym_id') == asym_id and r.get('ptnr1_label_asym_id') == chain) or \
                   (r.get('ptnr1_label_asym_id') == asym_id and r.get('ptnr2_label_asym_id') == chain):
                    partner_key = 'ptnr1' if r.get('ptnr2_label_asym_id') == asym_id else 'ptnr2'
                    comp_id_key = f'{partner_key}_label_comp_id'
                    auth_seq_id_key = f'{partner_key}_auth_seq_id'

                    # Verify if the residue is a modified residue
                    if 'pdbx_role' in r and r['pdbx_role'][0] == '?':
                        modified_res_message = check_modified_residues(
                            cif_data, chain, r.get(comp_id_key), r.get(auth_seq_id_key))
                    else:
                        modified_res_message = r.get('pdbx_role', 'Unknown role')
                    

                    covalent = f"{r.get(comp_id_key, 'UnknownComp')}{r.get(auth_seq_id_key, 'UnknownSeq')} {chain}"
                    return covalent, modified_res_message

    return covalent, modified_res_message



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
    

    
def search_comp_type(cif_data, id):
    """
    For standard polymer components, the type of the monomer.
    """

    if '_chem_comp' in cif_data:
        for r in cif_data._chem_comp.search('id',id).values():
            comp_type = r['type']    
        
        return comp_type   
     


def ligand_info(cif_data, entity_id, asym_id, entity_type, pdbx_description):
    """
    Search peptide ligand info.
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


def extract_sequences(fasta_file):
    """
    Extract sequences from a FASTA file as a dictionary {fasta_id: sequence}.
    """
    sequences_dict = {}
    
    if fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences_dict[str(record.id)] = str(record.seq)
    
    return sequences_dict



def search_for_mutation(sequences_dict, one_letter_seq):
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
    
    best_identity = 0
    best_alignment = None
    mismatches = None
    mismatch_location = None
    gaps = None
    best_reference_seq = None
    
    for fasta_id, reference_seq in sequences_dict.items():
        # Perform pairwise alignment
        alignments = aligner.align(reference_seq, one_letter_seq)
        alignment = alignments[0]  # We select the first alignment. It could be more than 1 alignment with same score
        identity = alignment.counts()[1] * 100 / len(alignment[0, :])
        
        if identity > best_identity:
            best_identity = identity
            best_alignment = alignment
            gaps = alignment.counts()[0]
            mismatches = alignment.counts()[2]
            best_reference_seq = reference_seq
            best_fasta_id = fasta_id  # Store the corresponding ID
            
            # Where is the mismatch
            n = 0
            mismatch_location = []

            for i, j in zip(alignment[0, :], alignment[1, :]):
                if i != '-':
                    n += 1
                    if i != j and j != '-':
                        mismatch_location.append(i + str(n) + j)
    
    identity = '{:.2f}'.format(best_identity)

    return mismatches, mismatch_location, identity, gaps, best_reference_seq, best_fasta_id


def read_coordinates_cif(file_path, entity_id):
    """
    Read coordinates and residue names.
    
    Returns:
        np.ndarray: coordinates (x, y, z) for each entity.
        list: residue names corresponding to each coordinate.
    """
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path)
    block = list(cif_obj.keys())[0]
    
    atom_site = cif_obj[block]["_atom_site"]
    
    # Filter by entity_id and skip hydrogen atoms
    mask = (np.array(atom_site["label_entity_id"]) == entity_id) & (np.array(atom_site["type_symbol"]) != "H")
    
    # Extract coordinates
    x = np.array(atom_site["Cartn_x"], dtype=float)[mask]
    y = np.array(atom_site["Cartn_y"], dtype=float)[mask]
    z = np.array(atom_site["Cartn_z"], dtype=float)[mask]
    
    # Check if any coordinates were found
    if len(x) == 0:
       
        return np.array([]), []
    
    # Extract atomic details
    atom_ids = np.array(atom_site["label_atom_id"])[mask]
    atom_types = np.array(atom_site["type_symbol"])[mask]
    residue_names = np.array(atom_site["label_comp_id"])[mask]
    residue_numbers = np.array(atom_site["label_seq_id"])[mask]
    occupancies = np.array(atom_site["occupancy"], dtype=float)[mask]

    # Combine all atomic details including occupancy
    atoms = list(zip(atom_ids, atom_types, residue_names, residue_numbers, occupancies))

    return np.column_stack((x, y, z)), atoms


def ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details):
    """
    Extract ligand information and create a dictionary with the ligand details.
    """
    
    if ligand_id not in ligands:
        ligands.append(ligand_id)

    if ligand_id not in ligand_details:
        ligand_details[ligand_id] = {}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        ligand_details[ligand_id]['type'] = ligand_type
        ligand_details[ligand_id]['function'] = ligand_function
        ligand_details[ligand_id]['name'] = entity_pdbx_description
        ligand_details[ligand_id]['entity'] = entity_id
        ligand_details[ligand_id]['covalent'] = ''
        ligand_details[ligand_id]['covalent_bond'] = ''

    if entity_type == 'non-polymer':
        comp_type = search_comp_type(cif_data, ligand_id)
        if comp_type: # If the _chem_comp.type exist, use the _chem_comp.type
            ligand_details[ligand_id]['type'] = comp_type
        else:
            ligand_details[ligand_id]['type'] = 'small-molecule'
     
    
    return ligand_details, ligands



def search_covalent_bonds(ligands, cif_data, ligand_details, protein_coords_map, ligand_coords_map, search_asym_id, chain, covalent_distance = 1.95):
    """
    Search for covalent bonds between ligands and proteins at the atomic level.
    
    Returns:
        Updated dictionary with covalent bonds.
    """
    
    
    for i in ligands:
        bond = ''
        covalent = 'No'
        asym_id = search_asym_id(cif_data, ligand_details[i]['entity'])

        for asym in asym_id.split(','):
            message = ''
            covalent_bond, message = check_covalent(cif_data, asym, chain)

            if covalent_bond:
                covalent = f'({message})' if message and message not in ["Unknown role", "covalent", "?"] else 'Yes'
                bond = bond + '; ' + covalent_bond if bond else covalent_bond
                break  # Exit if a covalent bond is found

        # # If no covalent bond was found, check atom-level distances
        if (covalent == 'No' and i in ligand_coords_map) or covalent == 'Yes':
            ligand_coordinates = ligand_coords_map[i]['coordinates']
            ligand_atoms = ligand_coords_map[i]['residues']  # Now contains atomic details
            best_distance = float('inf')
            best_bond = ''
            found_covalent = False

            for asym, protein_data in protein_coords_map.items():
                protein_coordinates = protein_data['coordinates']
                protein_atoms = protein_data['residues']  # Includes occupancy

                if protein_coordinates is None or not isinstance(protein_coordinates, np.ndarray) or protein_coordinates.size == 0:
                    continue

                distances = np.linalg.norm(
                    protein_coordinates[:, np.newaxis, :] - ligand_coordinates[np.newaxis, :, :], axis=2
                )
                is_close = distances < covalent_distance
                pairs = np.argwhere(is_close)

                if pairs.size > 0:
                    for pair in pairs:
                        p_idx, l_idx = pair
                        protein_atom = protein_atoms[p_idx]
                        ligand_atom = ligand_atoms[l_idx]
                        distance = distances[p_idx, l_idx]

                        # Comprobar occupancy
                        try:
                            protein_occupancy = float(protein_atom[4])
                            ligand_occupancy = float(ligand_atom[4])
                        except (IndexError, ValueError):
                            protein_occupancy = ligand_occupancy = 0.0

                        if protein_occupancy == 1.0 and ligand_occupancy == 1.0 and distance < best_distance:
                            best_distance = distance
                            found_covalent = True
                            best_bond = (
                                f'{protein_atom[0]} ({protein_atom[2]}{protein_atom[3]}) (Occup.: {protein_occupancy:.2f}) - '
                                f'{ligand_atom[0]} ({ligand_atom[2]}{ligand_atom[3]}) (Occup.: {ligand_occupancy:.2f}) - '
                                f'({distance:.2f} Å)'
                            )

            if found_covalent:
                covalent = 'Yes'
                bond = best_bond
            else:
                covalent = 'No'
                bond = ''

        ligand_details[i]['covalent'] = covalent
        ligand_details[i]['covalent_bond'] = bond
       
    return ligand_details


def process_cif_file(file_path, mutation, blacklist, sequences_dict, res_threshold, covalent_distance):
    """ 
    Process the cif file analyzing specific data as: ID, chain, number of residues, if it is a complex, 
    if it is peptide-like complex, the type of the bond and the name of the ligand
    """
     
    # Function to process a CIF file and extract the required information
    cfr = CifFileReader()
    cif_obj = cfr.read(file_path, output = 'cif_wrapper') 
    cif_data = list(cif_obj.values())[0]

    pdb_id = file_path.split('/')[-1].split('.')[0]
    title = cif_data['_struct']['title'][0]
    crystallization_method = cif_data['_exptl']['method'][0]
        
    # # ===================
    # # No NMR structures are analyzed
    # # ===================
    if "SOLUTION NMR" in crystallization_method or "SOLID-STATE NMR" in crystallization_method:
        print(f"Warning. NMR structure {pdb_id} is not analyzed.")
        return None
    
    # ===================
    # Initialize variables
    # ===================
    ligands = []
    ligand_names = []
    ligand_types = []
    ligand_functions = []
    ligand_covalents = []
    ligand_covalents_bond = []
    ligand_details = {}

    protein_description = []
    chain = '' 
    num_res = []  
    other_entities = []

    mismatches_list = []
    mismatch_location_list = []
    identity_list = [] 
    gaps_list = []
    best_fasta_id = None

    blacklist_id = []

    ligand_coords_map = {}
    protein_coords_map = {}
    protein_coordinates = None
    protein_residues = None

    # ===================
    # Search for entities
    # ===================
    for i in list(zip(cif_data['_entity']['id'], cif_data['_entity']['type'], cif_data['_entity']['src_method'], cif_data['_entity']['pdbx_description'])):
        
        entity_id = i[0]
        entity_type = i[1]
        entity_pdbx_description = i[3]
        asym_id = search_asym_id(cif_data, entity_id)

        # ===================
        # TYPE 1: POLYMER
        # ===================
        if entity_type == 'polymer': # if _entity.type is 'polymer' we deduce it is the protein or a peptide ligand  
            
            if check_polypeptide(cif_data, entity_id) == 'Yes': # 'polypeptide(L)' polymers could contain the protein or peptide ligands

                seq = search_peptide_seq_one_letter_code_can(cif_data, entity_id)
                if not seq:
                    print(f"Warning. No sequence found for entity_id {entity_id}")
                    continue   

                if len(seq) < res_threshold: # Small peptides are considered ligands
                    ligand_id, ligand_function, ligand_type = ligand_info(cif_data, entity_id, asym_id[0], entity_type, entity_pdbx_description)
                    ligand_details, ligands = ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details)
                    ligand_coordinates, ligand_residues = read_coordinates_cif(file_path, entity_id)
                    ligand_coords_map[ligand_id] = {'coordinates': ligand_coordinates,'residues': ligand_residues} 
                    
                # Proteins
                else:
                    chain = asym_id if chain == '' else chain + ',' + asym_id
                    protein_coordinates, protein_residues = read_coordinates_cif(file_path, entity_id)
                    protein_coords_map[asym_id] = {'coordinates': protein_coordinates,'residues': protein_residues} 
                    protein_description.append(entity_pdbx_description)
                    num_res.append(str(len(seq)))
                    
                    mismatches, mismatch_location, identity, gaps, best_reference_seq, best_fasta_id = (None, None, None, None, None, None)
        
                    if mutation == True:
                        mismatches, mismatch_location, identity, gaps, best_reference_seq, best_fasta_id = search_for_mutation(sequences_dict, seq)

                    mismatches_list.append(mismatches if mismatches is not None else "")
                    mismatch_location_list.append(mismatch_location if mismatch_location is not None else [])
                    identity_list.append(identity if identity is not None else "")
                    gaps_list.append(gaps if gaps is not None else "")

            else: # 'polypeptide(D)' or DNA/RNA
                ligand_id, ligand_function, ligand_type = ligand_info(cif_data, entity_id, asym_id[0], entity_type, entity_pdbx_description)
                ligand_details, ligands = ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details)
                ligand_coordinates, ligand_residues = read_coordinates_cif(file_path, entity_id)
                ligand_coords_map[ligand_id] = {'coordinates': ligand_coordinates,'residues': ligand_residues} 
            
        # ===================       
        # TYPE 2: NON-POLYMER
        # ===================
        elif entity_type == 'non-polymer': # non-polymer could be ligands or small-molecules from the medium
            
            ligand_id = search_ligand_code(cif_data, entity_id) # From the entity_id, search the ligand_code or sequence (for peptides without a code)
            
            if ligand_id not in blacklist: # Only molecules not present in the blacklist file are considered as ligands
                ligand_type = 'non-polymer'
                ligand_function = ''
                ligand_details, ligands = ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details)
                ligand_coordinates, ligand_residues = read_coordinates_cif(file_path, entity_id)  
                ligand_coords_map[ligand_id] = {'coordinates': ligand_coordinates,'residues': ligand_residues}      
                
            else:
                blacklist_id.append(ligand_id)
        
        # ===================
        # TYPE 3: BRANCHED 
        # ===================
        if entity_type == 'branched':
            pep_code = search_peptide_code(cif_data, asym_id[0])
            
            if pep_code:
                ligand_id, ligand_function, ligand_type = ligand_info(cif_data,entity_id,asym_id[0],entity_type,entity_pdbx_description)
                ligands.append(ligand_id)
                ligand_details, ligands = ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details)
                ligand_coordinates, ligand_residues = read_coordinates_cif(file_path, entity_id)
                ligand_coords_map[ligand_id] = {'coordinates': ligand_coordinates,'residues': ligand_residues}
            
            else:
                ligand_id, ligand_type = search_branched(cif_data, entity_id) 
                ligands.append(ligand_id)
                ligand_function = ''
                ligand_details, ligands = ligand_dictionary(cif_data, entity_type, ligand_id, ligand_function, ligand_type, entity_id, entity_pdbx_description, ligands, ligand_details)
                ligand_coordinates, ligand_residues = read_coordinates_cif(file_path, entity_id)
                ligand_coords_map[ligand_id] = {'coordinates': ligand_coordinates,'residues': ligand_residues}
    
    # Search for covalent bonds
    search_covalent_bonds(ligands, cif_data, ligand_details, protein_coords_map, ligand_coords_map, search_asym_id, chain, covalent_distance)

    # Ligand details dictionary
    for i in ligand_details.keys():
        ligand_names.append(ligand_details[i]['name'])
        ligand_types.append(ligand_details[i]['type'])
        ligand_functions.append(ligand_details[i]['function'])
        ligand_covalents.append(ligand_details[i]['covalent'])
        ligand_covalents_bond.append(ligand_details[i]['covalent_bond'])

    ligand_names = '\n'.join(ligand_names)
    ligand_types = '\n'.join(ligand_types)
    ligand_functions = '\n'.join(ligand_functions)
    ligand_covalents = '\n'.join(ligand_covalents)
    ligand_covalents_bond = '\n'.join(ligand_covalents_bond) # igual aqui
    ligands = '\n'.join(ligands)
    complex = 'Yes' if ligands else 'No'

    # Convert lists to strings for the output
    protein_description = '\n'.join(protein_description)
    nchain = chain.count(',') + 1
    num_res = ', '.join(num_res)
    blacklist_id = '\n'.join(blacklist_id)
    other_entities = '\n'.join(other_entities)
    fasta_id = best_fasta_id if mutation else None
    mismatch_location = [", ".join(map(str, subunit)) for subunit in mismatch_location_list] if mutation else []
    mismatch_location = '\n'.join(mismatch_location)
    mismatches = '\n'.join([str(m) for m in mismatches_list if m is not None])
    identity = '\n'.join([str(id) for id in identity_list if id is not None])
    gaps = '\n'.join([str(g) for g in gaps_list if g is not None])
    

    return {
        "PDB_ID": pdb_id,
        "Title": title,
        "Protein_description": protein_description,
        "NChain": nchain,
        "Chain_ID": chain,
        "Num_Res": num_res,
        "Complex": complex,
        "Discarded_Ligands": blacklist_id,
        "Ligand": ligands,
        "Ligand_names": ligand_names,
        "Ligand_types": ligand_types,
        "Ligand_functions": ligand_functions,
        "Covalent_Bond": ligand_covalents,
        "Bond": ligand_covalents_bond,
        "FASTA_ID": fasta_id,
        "Mutation": mismatches,
        "Mutation_Location": mismatch_location,
        "Identity": identity,
        "Gaps": gaps
    }


def _process_files(file_list, mutation, blacklist, sequences_dict, res_threshold, blacklist_dict, fields_to_include, covalent_distance):
    """
    Processes a single .cif file and returns extracted data.
    """
    global valid_count
    global processed_count
    data, data_ligands = [], []
    
    for file_path in file_list:
        try:
            data_from_file = process_cif_file(file_path, mutation, blacklist, sequences_dict, res_threshold, covalent_distance)
            with processed_count.get_lock():
                processed_count.value += 1
            if data_from_file is None:
                continue

            ligands = data_from_file.get("Ligand", "").split('\n')
            ligand_names_list = data_from_file.get("Ligand_names", "").split('\n')
            ligand_types_list = data_from_file.get("Ligand_types", "").split('\n')
            covalent_bond_list = data_from_file.get("Covalent_Bond", "").split('\n')
            ligand_covalents_bond = data_from_file.get("Bond", "").split('\n')
            discarded_ligands = data_from_file.get("Discarded_Ligands", "").split('\n')
            
            max_length = max(map(len, [ligands, ligand_names_list, ligand_types_list, covalent_bond_list, ligand_covalents_bond, discarded_ligands]))
            ligand_data = [
                {
                    **{field: data_from_file.get(field, "") for field in fields_to_include},
                    "Ligand": ligands[i].strip() if i < len(ligands) else "",
                    "Ligand_names": ligand_names_list[i].strip() if i < len(ligand_names_list) else "",
                    "Ligand_types": ligand_types_list[i].strip() if i < len(ligand_types_list) else "",
                    "Covalent_Bond": covalent_bond_list[i].strip() if i < len(covalent_bond_list) else "",
                    "Bond": ligand_covalents_bond[i].strip() if i < len(ligand_covalents_bond) else "",
                }
                for i in range(max_length)
            ]
            
            ligand_data += [
                {
                    **{field: data_from_file.get(field, "") for field in fields_to_include},
                    "Ligand": discarded_ligands[i].strip(),
                    "Ligand_names": blacklist_dict.get(discarded_ligands[i].strip(), ""),
                    "Ligand_types": "Discarded",
                    "Covalent_Bond": "",
                    "Bond": "",
                }
                for i in range(len(discarded_ligands)) if discarded_ligands[i].strip()
            ]
            
            data.append(data_from_file)
            data_ligands.extend(ligand_data)
            
            with valid_count.get_lock():
                valid_count.value += 1 
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    return data, data_ligands



def write_output(directory_path, out_file, out_file_ligands, blacklist_dict, mutation, blacklist, sequences_dict, res_threshold, covalent_distance=1.95, cpu_use=1.0):
    """
    Processes .cif files in parallel and generates two optimized CSV files, with system monitoring.
    """
    
    def monitor_system(total_files):
        """
        Monitors CPU, RAM, and file processing progress in the background.
        """
        start_time = time.time()
        
        while not stop_event.is_set():
            elapsed_time = time.time() - start_time
            
            mem_info = psutil.virtual_memory()
            ram_usage_mb = mem_info.used / 1024**2   # Convert to MB
            ram_usage_percent = mem_info.percent
            cpu_usage = psutil.cpu_percent(interval=0.5, percpu=True)  # Per-core CPU usage
            completed = processed_count.value

            # First line: General stats
            info_line = (
                f"Time: {elapsed_time:.2f}s | "
                f"RAM Usage: {ram_usage_mb:.2f} MB ({ram_usage_percent:.2f}%) | "
                f"CPU: {sum(cpu_usage) / len(cpu_usage):.2f}% | "
                f"Files: {completed}/{total_files}"
            )
            
            # Second line: Per-core CPU usage
            cpu_line = "CPU Usage per core: " + " | ".join(f"{i + 1}: {usage:.1f}%" for i, usage in enumerate(cpu_usage))

            # Clear previous output and rewrite lines
            global impression
            if impression == True:
                print("\033[2J\033[H", end="")  # ANSI escape codes to clear screen
            else: impression = True
            print(info_line)
            print(cpu_line)

            time.sleep(1)

    stop_event = threading.Event()  
    cif_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith(('.cif', '.cif.gz'))]
    total_files = len(cif_files)
    if total_files == 0:
        print("No .cif or .cif.gz files found in the directory.")
        return
    
    monitor_thread = threading.Thread(target=monitor_system, args=(total_files,), daemon=True)  
    monitor_thread.start()  

    # ============================
    # CONFIGURATION
    # ============================
    
    fields_to_include = ["PDB_ID", "Ligand", "Ligand_names", "Ligand_types", "Ligand_functions", "Covalent_Bond", "Bond"]
    num_cpus = min(total_files, max(1,int(multiprocessing.cpu_count()*min(1,cpu_use))))
    chunk_size = max(1, total_files // num_cpus)
    file_chunks = [cif_files[i:i + chunk_size] for i in range(0, total_files, chunk_size)]

    data, data_ligands = [], []

    # ============================
    # PARALLEL PROCESSING
    # ============================

    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = {
            executor.submit(_process_files, chunk, mutation, blacklist, sequences_dict, res_threshold, blacklist_dict, fields_to_include, covalent_distance): chunk 
            for chunk in file_chunks
        }
        
        for future in as_completed(futures):
            try:
                chunk_data, chunk_ligands = future.result()
                data.extend(chunk_data)
                data_ligands.extend(chunk_ligands)
            except Exception as e:
                print(f"Error processing chunk: {e}")

    # ============================
    # OPTIMIZED FILE SAVING
    # ============================

    pd.DataFrame(data).to_csv(out_file, index=False)

    df_ligand = pd.DataFrame(data_ligands)
    df_ligand['Ligand'] = df_ligand['Ligand'].str.strip()
    df_ligand = df_ligand[df_ligand['Ligand'] != '']
    
    df_ligand.columns = ['ID', 'Molecule', 'Name', 'Type', 'Function', 'Covalent', 'Bond']
    df_ligand.to_csv(out_file_ligands, index=False)

    # ============================
    # STOP SYSTEM MONITORING
    # ============================

    stop_event.set()  
    monitor_thread.join()  
    print(f"Files successfully analyzed: {valid_count.value}/{total_files}")

    return data, data_ligands


def mutation_classification(directory_path, information_df, output_path):
    """
    Classification based on the information contained in a dataframe. 
    Downloading PDB files into different folders depending on whether there is a mutation.
    """
    mutated_list = []
    no_mutated_list = []
    
    df = pd.read_csv(information_df, dtype={'Mutation': str})
    
    for _, row in df.iterrows():
        if row['Mutation'] == '0':
            no_mutated_list.append(row['PDB_ID'])
        else:
            mutated_list.append(row['PDB_ID'])
    
    # Generate timestamped folders if needed
    current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    paths = {
        "mutated": os.path.join(output_path, "Mutated"),
        "non_mutated": os.path.join(output_path, "Non-Mutated")
    }
    
    for key in paths:
        if os.path.exists(paths[key]):
            paths[key] = f"{paths[key]}_{current_datetime}"
        os.makedirs(paths[key], exist_ok=True)
    
    def copy_pdb(pdb_id, dest_path):
        for ext in ['.cif', '.cif.gz']:
            src = os.path.join(directory_path, pdb_id + ext)
            dest = os.path.join(dest_path, pdb_id + ext)
            if os.path.exists(src):
                shutil.copy(src, dest)
    
    for pdb_id in mutated_list:
        copy_pdb(pdb_id, paths["mutated"])
    
    for pdb_id in no_mutated_list:
        copy_pdb(pdb_id, paths["non_mutated"])
    
    return no_mutated_list, paths["non_mutated"]


def bond_classification(directory_path, information_df, no_mutated_list, output_path, mutation):
    """
    Classifies between free enzymes and complexes.
    For each complex, classifies based on their bonds: covalent or non-covalent.
    """
    
    free_enzyme_list = []
    covalent_list = []
    non_covalent_list = []
    
    df = pd.read_csv(information_df)
    
    for _, row in df.iterrows():
        if row['PDB_ID'] in no_mutated_list:
            if row['Complex'] == 'No':
                free_enzyme_list.append(row['PDB_ID'])
            else:
                if isinstance(row['Covalent_Bond'], str) and 'Yes' in row['Covalent_Bond']:
                    covalent_list.append(row['PDB_ID'])
                else:
                    non_covalent_list.append(row['PDB_ID'])
    
    origin_path = output_path if mutation else directory_path
    
    paths = {
        "free": os.path.join(output_path, 'APO'),
        "covalent": os.path.join(output_path, 'Covalent'),
        "non_covalent": os.path.join(output_path, 'Non-Covalent')
    }
    
    # Create directories if they don't exist
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
    
    def move_or_copy(pdb_id, dest_path):
        for ext in ['.cif', '.cif.gz']:
            src = os.path.join(origin_path, pdb_id + ext)
            dest = os.path.join(dest_path, pdb_id + ext)
            if os.path.exists(src):
                if mutation:
                    shutil.move(src, dest)
                else:
                    shutil.copy(src, dest)
    
    for pdb_id in free_enzyme_list:
        move_or_copy(pdb_id, paths["free"])
    
    for pdb_id in covalent_list:
        move_or_copy(pdb_id, paths["covalent"])
    
    for pdb_id in non_covalent_list:
        move_or_copy(pdb_id, paths["non_covalent"])
    
    return
