import pandas as pd
import os
from paths_and_settings import *
from pdb_query import get_pdbs_from_unicode, get_pdb_fasta, check_if_pdb_xray
from utils import read_fasta
from urllib.error import HTTPError


def choose_primary_pdb_for_chembl(csv_path):
    main_table = pd.read_csv(csv_path, index_col=0)
    main_table['main_PDB_structure'] = ''
    for index, row in main_table.iterrows():
        pdbs = row['PDB_entry'].split()
        fasta = [get_pdb_fasta(i) for i in pdbs]
        best = ''
        for seq in fasta:
            if len(seq) > len(best):
                best = seq
        if best != '':
            best_pdb = pdbs[fasta.index(best)]
            main_table.at[index, 'main_PDB_structure'] = best_pdb
    main_table.to_csv(csv_path)


def fetch_fastas_for_chembl(csv_path, output):
    main_table = pd.read_csv(csv_path, index_col=0)
    print(f'Downloading fastas for {len(main_table)} targets')
    with open(output, 'w') as handle:
        for index, row in main_table.iterrows():
            chembl_name = row['ChEMBL ID']
            pdb_name = row['main_PDB_structure']
            fasta = get_pdb_fasta(row['main_PDB_structure'])
            handle.write(f'>{chembl_name}-{pdb_name}\n{fasta}\n')


def fetch_fastas_for_DEKOIS(dekois_dir_path=DEKOIS_PATH):
    print('Starting fetching FASTA\' s for DEKOIS')
    dekois_ligands, _ = DEKOIS_uniID_from_folder(dekois_dir_path)
    dekois_ligands_uni_ids = set()
    for i in dekois_ligands:
        try:
            for k in dekois_ligands[i]['Uniprot_ID'].split():
                dekois_ligands_uni_ids.add(k)
        except KeyError:
            print(f'No Uniprot ID for {i}')
    with open('fastas_from_dekois.txt', 'w') as handle:
        for uniprot_id in dekois_ligands_uni_ids:
            pdb_ids = get_pdbs_from_unicode(uniprot_id)
            pdbs_dict = dict()
            for pdb in pdb_ids:
                try:
                    if check_if_pdb_xray(pdb):
                        pdb_fasta = get_pdb_fasta(pdb)
                        pdbs_dict[uniprot_id + '-' + pdb] = pdb_fasta.strip('"')
                except HTTPError:
                    print(f'Failed fetching FASTA sequence for: {pdb}')
            pdbs_dict2 = dict()
            for key, value in pdbs_dict.items():
                if value not in list(pdbs_dict2.values()):
                    pdbs_dict2[key] = value
            for key, value in pdbs_dict2.items():
                handle.write(f'>{key}\n{value}\n')


def fetch_fastas_for_DUDE(dude_dir_path=DUDE_PATH):
    print(f'Starting fetching FASTA\'s for DUDE')
    dude_uni_ids = DUDE_uniID_from_folder(dude_dir_path)
    with open('fastas_from_dude.txt', 'w') as handle:
        for uniprot_id in dude_uni_ids:
            pdb_ids = get_pdbs_from_unicode(uniprot_id)
            pdbs_dict = dict()
            for pdb in pdb_ids:
                try:
                    if check_if_pdb_xray(pdb):
                        pdb_fasta = get_pdb_fasta(pdb)
                        pdbs_dict[uniprot_id + '-' + pdb] = pdb_fasta.strip('"')
                except HTTPError:
                    print(f'Failed fetching FASTA sequence for: {pdb}')
            pdbs_dict2 = dict()
            for key, value in pdbs_dict.items():
                if value not in list(pdbs_dict2.values()):
                    pdbs_dict2[key] = value
            for key, value in pdbs_dict2.items():
                handle.write(f'>{key}\n{value}\n')


def DUDE_uniID_from_folder(dude_folder=DUDE_PATH):
    dude_uni_ids = set()
    for dirpath, dirnames, filenames in os.walk(dude_folder):
        for item in filenames:  # loop through items in dir
            if item.endswith('uniprot.txt'):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of file
                with open(file_name, 'r') as uni_file:
                    for line in uni_file:
                        dude_uni_ids.add(line.strip())
    return dude_uni_ids


def DEKOIS_uniID_from_folder(dir_name, extension='.sdf',
                             search_terms=('Name', 'Uniprot_ID', 'IC50_in_nM', 'Ki_in_nM', 'Kd_in_nM')):
    decoy_dir = os.path.join(dir_name, 'decoys')
    ligand_dir = os.path.join(dir_name, 'ligands')
    uniprot_ligands_ids = dict()
    uniprot_decoys_ids = dict()
    for dirpath, dirnames, filenames in os.walk(ligand_dir):
        for item in filenames:  # loop through items in dir
            if item.endswith(extension):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of files
                with open(file_name, 'r') as sdf_file:
                    for line in sdf_file:
                        for search_term in search_terms:
                            if search_term in line:
                                if search_term == 'Name':
                                    next_line = sdf_file.readline().strip()
                                    uniprot_ligands_ids[next_line] = dict()
                                else:
                                    uniprot_ligands_ids[next_line][search_term] = sdf_file.readline().strip()
    for dirpath, dirnames, filenames in os.walk(decoy_dir):
        for item in filenames:  # loop through items in dir
            if item.endswith(extension):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of files
                with open(file_name, 'r') as sdf_file:
                    for line in sdf_file:
                        for search_term in search_terms:
                            if search_term in line:
                                if search_term == 'Name':
                                    next_line = sdf_file.readline().strip()
                                    uniprot_decoys_ids[next_line] = dict()
                                else:
                                    uniprot_decoys_ids[next_line][search_term] = sdf_file.readline().strip()

    with open(os.path.join(BLAST_MAIN_FOLDER, 'dekois_ligands_uniprot.txt'), 'w') as handle:
        for target in uniprot_ligands_ids:
            handle.write(target + '\n')
            for item in uniprot_ligands_ids[target]:
                handle.write(item + '\t' + uniprot_ligands_ids[target][item] + '\n')

    with open(os.path.join(BLAST_MAIN_FOLDER, 'dekois_decoys_uniprot.txt'), 'w') as handle:
        for target in uniprot_decoys_ids:
            handle.write(target + '\n')
            for item in uniprot_decoys_ids[target]:
                handle.write(item + '\t' + uniprot_ligands_ids[target][item] + '\n')

    return uniprot_ligands_ids, uniprot_decoys_ids


def read_pdbbind_ids(file):
    pdb_ids = []
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '#':
                pdb_ids.append(line.split(' ')[0].upper())
    return pdb_ids


def load_fasta_sequences(fasta_file1, fasta_file2, fasta1_to_fasta2_blast):
    fasta_sequences = {''.join(filename.split('.')[:-1]): [[headerStr, seq] for headerStr, seq in read_fasta(filename)]
                       for filename in [fasta_file1, fasta_file2]}
    blast_results = load_blast_results(fasta1_to_fasta2_blast, [0, 1, 2])
    fasta_sources = list(fasta_sequences.keys())
    output_name = '-'.join([fasta_sources[0], fasta_sources[1]])
    cols = [v[0] for v in fasta_sequences[fasta_sources[0]]]
    rows = [v[0] for v in fasta_sequences[fasta_sources[1]]]
    values = pd.DataFrame(0.0, index=rows, columns=cols)
    for val in blast_results:
        values.at[val[1], val[0]] = val[2]
    values.to_csv(output_name + '.csv')


def load_blast_results(blast_results, usecols):
    df = pd.read_csv(blast_results, header=None, usecols=usecols)
    values = df.values.tolist()
    return values


def make_blast_csv(master_path, blasts=None):
    main_table = pd.read_csv(master_path, index_col=0)
    output_df = main_table.loc[:, ['ChEMBL ID', 'main_PDB_structure', 'Active_compounds', 'Inactive_compounds']]
    output_file = os.path.join(BLAST_MAIN_FOLDER, 'chembl_blast_results.csv')
    for db in blasts:
        #chembl_Id, hit_ID, identity%, evalue
        blast_results = load_blast_results(db, [0, 1, 2, 6, 7, 10])
        db_name = str(os.path.split(db)[-1]).split('_')[0].split('-')[1]
        output_df[f'identity%_{db_name}'] = ''
        output_df[f'evalue_{db_name}'] = ''
        output_df[f'target_name_{db_name}'] = ''
        output_df[f'query_alignment_length_{db_name}'] = ''
        output_df[f'total_query_length_{db_name}'] = ''
        output_df[f'alignment_to_total_ratio_{db_name}'] = ''
        name = ''
        for chembl_id, db_target, identity, q_start, q_end, evalue in blast_results:
            if name not in chembl_id or name == '':
                query_pdb_smile_len = len(get_pdb_fasta(chembl_id.split('-')[-1]))
                query_alignment_length = q_end - q_start
                name = chembl_id.split('-')[0]
                output_df.at[name, f'identity%_{db_name}'] = identity
                output_df.at[name, f'evalue_{db_name}'] = evalue
                output_df.at[name, f'target_name_{db_name}'] = db_target
                output_df.at[name, f'query_alignment_length_{db_name}'] = query_alignment_length
                output_df.at[name, f'total_query_length_{db_name}'] = query_pdb_smile_len
                output_df.at[name, f'alignment_to_total_ratio_{db_name}'] = query_alignment_length/query_pdb_smile_len
    output_df.to_csv(output_file)
