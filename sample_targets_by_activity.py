import pandas as pd
import pypdb
import pickle
from collections import Counter
from rdkit import DataStructs
import os
import numpy as np
from utils import make_fingerprint, log
from pdb_query import get_pdb_fasta
from paths_and_settings import *


def get_ligands(pdbs, suffix=""):
    results = dict()
    pdb_to_ligands = dict()
    s = 1
    for target in pdbs:
        if s % 100 == 0:
            log(f'Currently checked {s} targets.')
        ligands = []
        smiles = dict()
        ex = pdbs[target].split()
        log(f'{len(ex)} PDB\'s found for target {target}, {s}/{len(pdbs)}')
        k = 0
        failes = 0
        for pdb_id in ex:
            while True:
                try:
                    pdb_to_ligands[pdb_id] = []
                    i_ligands = pypdb.get_ligands(pdb_id)
                    if i_ligands['ligandInfo'] is not None:
                        ligs = i_ligands['ligandInfo'][
                            'ligand']  # l['ligandInfo']['ligand'] --> '@structureId', '@chemicalID', 'smiles'
                        for lig in ligs:
                            if type(lig) == dict:
                                ligands.append(str(lig['@chemicalID']))
                                pdb_to_ligands[pdb_id].append(str(lig['@chemicalID']))
                                smiles[lig['@chemicalID']] = lig['smiles']
                    failes = 0
                except:
                    failes += 1
                    if failes > 200:
                        log(f'Can\'t download ligand {pdb_id}, retried {failes} times.')
                    continue
                break
            k += 1
        results[target] = [ligands, smiles]
        s += 1
    pickle.dump(results, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"ligands{suffix}.pkl"), "wb"))
    pickle.dump(pdb_to_ligands, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"pdb_to_ligands{suffix}.pkl"), "wb"))
    return results, pdb_to_ligands


def count_ligands(results):
    targets = list(results.keys())
    res = []
    for t in targets:
        ligs = results[t][0]
        res = res + ligs
    res = dict(Counter(res))
    return res


def make_fingerprints_base(targets, index='', path=STANDARDS_FOLDER, prefix='compounds_'):
    files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and prefix in i]
    standard_csvs = [pd.read_csv(os.path.join(path, i), low_memory=False)[['ID_compound', 'ID_target_CHEMBL', 'Target_SMILES']].dropna()
                     for i in files]
    pdbs = list(targets['ChEMBL ID'])
    smile_dict = dict()
    fails = []
    for target in pdbs:
        smile_dict[target] = []
        p = set()
        s2 = dict()
        for standard in standard_csvs:
            target_standard_smiles = set(standard[standard['ID_target_CHEMBL'] == target]['Target_SMILES'])
            p.update(target_standard_smiles)
            s3 = standard[standard['ID_target_CHEMBL'] == target][['ID_compound', 'Target_SMILES']]
            for index2, row in s3.iterrows():
                s2[row['Target_SMILES']] = row['ID_compound']
        p = list(p)
        k_inner = []
        for sm in p:
            try:
                chem_smiles = make_fingerprint(sm, fp='fp2')
                k_inner.append((chem_smiles, sm, s2[sm]))
            except:
                print('Failed to create Chem SMILES from: ', sm)
                fails.append((target, sm))
        smile_dict[target] = k_inner

    name_fp = os.path.join(COMPOUND_SAMPLING_FOLDER, f"fingerprints_chembl{index}.pkl")
    name_failes = os.path.join(COMPOUND_SAMPLING_FOLDER, f"fingerprints_chembl_fails{index}.pkl")
    pickle.dump(smile_dict, open(name_fp, "wb"))
    pickle.dump(fails, open(name_failes, "wb"))
    print(f'Saved fingerptints to {name_fp}')
    return smile_dict, fails


def pdb_ligands_vs_chembl_ligands(pdbs, fingerprints, threshold=0.95, suffix=""):
    targets = list(pdbs.keys())
    output = dict()
    fails = []
    for target in targets:
        output[target] = dict()
        target_smiles = list(pdbs[target][1].items())  # smiles ligandów z PDB
        target_chembl = fingerprints[target]
        for ligand_tuple in target_smiles:
            try:
                ligand_id = ligand_tuple[0]  # id ligandu
                ligand_smile = ligand_tuple[1]  # smiles ligandu
                output[target][ligand_id] = []
                l2 = make_fingerprint(ligand_smile, fp='fp2')
                for chembl_smile in target_chembl:
                    coef = DataStructs.FingerprintSimilarity(l2, chembl_smile[0])
                    if coef > threshold:
                        output[target][ligand_id].append((coef, chembl_smile[2]))
            except:
                fails.append((target, ligand_tuple))
        for i in list(output[target].keys()):
            if len(output[target][i]) == 0:
                output[target].pop(i)
    for i in list(output.keys()):
        if len(output[i]) == 0:
            output.pop(i, None)
    fails = [(i[0], i[1][0], i[1][1]) for i in fails]

    pickle.dump(output, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"pdb_vs_chembl_filtered_ligands_tc{threshold}-{suffix}.pkl"), "wb"))
    pickle.dump(fails, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"/fail_smiles_pdb_tc{threshold}-{suffix}.pkl"), "wb"))
    return output, fails


def filter_ligands(all_ligands, threshold=100):  # this function here is probably the biggest weak point of whole pipeline
    to_remove = set(list({k: v for k, v in count_ligands(all_ligands).items() if v > threshold}.keys()))
    for target in list(all_ligands.keys()):
        keys = set(all_ligands[target][0])
        intersection = keys.intersection(to_remove)
        for lig_name in intersection:
            all_ligands[target][0].remove(lig_name)
            all_ligands[target][1].pop(lig_name, None)
        if len(all_ligands[target][0]) == 0:
            all_ligands.pop(target, None)
    return all_ligands


def evaluate(suffix, tcs_list=(0.95,)):
    '''
    Make a matrix with comparison between runs on different parameters (how many proteins got through filtering).
    '''
    matrix = np.zeros((len(tcs_list), 1))
    for tc in range(len(tcs_list)):
        name = os.path.join(COMPOUND_SAMPLING_FOLDER, f'pdb_vs_chembl_filtered_ligands_tc{tcs_list[tc]}-{suffix}.pkl')
        with open(name, 'rb') as handle:
            filtered_ligands = pickle.load(handle)
        pickle.dump(filtered_ligands, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"targets_before_blast_tc{tcs_list[tc]}-{suffix}.pkl"), "wb"))
        with open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"targets_before_blast_tc{tcs_list[tc]}-{suffix}.txt"), 'w') as save_file:
            for i in filtered_ligands:
                save_file.write(i + '\n')
        matrix[tc, 0] = len(filtered_ligands)
    pickle.dump(matrix, open(os.path.join(COMPOUND_SAMPLING_FOLDER, f"runs_comparison_matrix{suffix}.pkl"), "wb"))
    pd.DataFrame(data=matrix, index=tcs_list).to_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, f'runs_comparison_matrix{suffix}.csv'))


def zip_results_with_pdb(results_final, pdb_to_ligands, master_table):
    output = dict()
    for target in results_final:
        output[target] = []  # (pdb_id, ligand, TC)
        try:
            possible_pdbs = list(master_table[master_table['ChEMBL ID'] == target]['PDB_entry'])[0].split()
        except:
            print('No PDB entrys for: ', target)
            continue
        target_ligands = list(results_final[target].keys())  # ligand
        for ligand in target_ligands:
            for pdb_id in possible_pdbs:
                try:
                    if ligand in pdb_to_ligands[pdb_id]:
                        output[target].append((pdb_id, ligand, results_final[target][ligand]))
                except:
                    pass
    return output


def update_master_table(master_table, results):
    master_table = master_table[master_table['ChEMBL ID'].isin(list(results.keys()))].copy(deep=True)
    master_table['PDBs_number'] = 0
    master_table['Ligand_number'] = 0
    for target in results:
        pdbs = set()
        ligs = set()
        for each in results[target]:
            pdbs.add(each[0])
            ligs.add(each[1])
        master_table.at[target, 'PDBs_number'] = len(pdbs)
        master_table.at[target, 'Ligand_number'] = len(ligs)
        new_pdbs = " ".join(list(pdbs))
        master_table.at[target, 'PDB_entry'] = new_pdbs
    return master_table


def choose_primary_pdb_for_chembl(main_table):
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
    main_table = main_table.sort_values('Active_compounds', ascending=False).drop_duplicates(subset='main_PDB_structure')
    return main_table


def make_chembls_smiles_files(path_to_fitered_targets, smiles_folder=STANDARDS_FOLDER, prefix='compounds', folder=CHEMBL_SMILES_FOLDER):
    """
    Make SMILES file for each chembl id from selected file in pickle format.
    :param path_to_fitered_targets: Pickled list of ids.
    :param smiles_folder: Where to find smiles csv's.
    :param prefix: Prefix of csv's with smiles.
    :param folder: Where to save smiles files.
    :return:
    """
    with open(path_to_fitered_targets, 'rb') as handle:
        filtered_targets = pickle.load(handle)
    files = [i for i in os.listdir(smiles_folder) if os.path.isfile(os.path.join(smiles_folder, i)) and prefix in i] # ZAŁADOWANE TARGET CHEMBL ID, ZNALEŹĆ ICH SMILES Z CHEMBL
    standard_csvs = [pd.read_csv(smiles_folder + '/' + i, low_memory=False).dropna() for i in files]
    smile_dict = dict()
    for chembl_id in filtered_targets:
        smile_dict[chembl_id] = dict()
        with open(f'{folder}/{chembl_id}_all.smi', 'w') as save_file:
            ligands = set()
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_target_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        smile = record[1]['Target_SMILES']
                        chembl_ligand = record[1]['ID_compound']
                        if chembl_ligand not in ligands:
                            save_file.write(smile + '\t' + chembl_ligand + '\n')
                            smile_dict[chembl_id][chembl_ligand] = smile
                            ligands.add(chembl_ligand)
                except KeyError:
                    pass
        with open(f'{folder}/{chembl_id}_active.smi', 'w') as save_file:
            ligands = set()
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_target_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        if record[1]['Activity'] == 'Active':
                            smile = record[1]['Target_SMILES']
                            chembl_ligand = record[1]['ID_compound']
                            if chembl_ligand not in ligands:
                                ligands.add(chembl_ligand)
                                save_file.write(smile + '\t' + chembl_ligand + '\n')
                except KeyError:
                    pass
        with open(f'{folder}/{chembl_id}_inactive.smi', 'w') as save_file:
            ligands = set()
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_target_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        if record[1]['Activity'] != 'Active':
                            smile = record[1]['Target_SMILES']
                            chembl_ligand = record[1]['ID_compound']
                            if chembl_ligand not in ligands:
                                ligands.add(chembl_ligand)
                            save_file.write(smile + '\t' + chembl_ligand + '\n')
                except KeyError:
                    pass
    return smile_dict


def sample_by_active_compounds(list_of_values):
    if not os.path.exists(COMPOUND_SAMPLING_FOLDER):
        os.makedirs(COMPOUND_SAMPLING_FOLDER)
    for val in list_of_values:
        targets = pd.read_csv(MAIN_CSV_NAME, sep=',', index_col=0)
        new = targets[pd.isna(targets['DEKOIS_ID']) & pd.isna(targets['DUDE_ID'])]
        new = new.sort_values(by=['Active_compounds'], ascending=False)
        new = new[new['Active_compounds'] >= val]
        new = new.drop(labels=['DEKOIS_ID', 'DEKOIS_actives',
       'DEKOIS_actives_threshold', 'DEKOIS_decoys', 'DUDE_ID',
       'Active_in_DUDE', 'Active_in_DUDE_threshold', 'Inactive_in_DUDE',
       'Inactive_in_DUDE_threshold'], axis=1)
        print(f'Targets without coverage for threshold {val}:', len(new))
        new.to_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, f'targets_after_filtering_{val}.csv'))


def sample_by_activity(activity_threshold, thresholds=(0.95,), ligand_threshold=100):
    sample_by_active_compounds(activity_threshold)
    for threshold in thresholds:
        for ix, file in enumerate([i for i in os.listdir(COMPOUND_SAMPLING_FOLDER) if 'targets_after_filtering_' in i]):
            index = file.split('.')[0].split("_")[-1]
            data = pd.read_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, file), index_col=0)
            pdbs = dict(data['PDB_entry'])
            try:
                with open(os.path.join(COMPOUND_SAMPLING_FOLDER, f'ligands{index}.pkl'), 'rb') as handle:
                    results = pickle.load(handle)
                with open(os.path.join(COMPOUND_SAMPLING_FOLDER, f'pdb_to_ligands{index}.pkl'), 'rb') as handle:
                    pdb_to_ligands = pickle.load(handle)
            except FileNotFoundError:
                results, pdb_to_ligands = get_ligands(pdbs, suffix=str(index))  # download ligands
                print(f'For index {index} - Number of targets with ligands: {len(results)}')
            try:
                with open(os.path.join(COMPOUND_SAMPLING_FOLDER, 'fingerprints_chembl.pkl'), 'rb') as handle:
                    smiles = pickle.load(handle)
                    print('Fingerprints found and loaded!')
            except FileNotFoundError:
                print('Preparing fingerprint base...')
                smiles, fails1 = make_fingerprints_base(targets=pd.read_csv(MAIN_CSV_NAME, index_col=0))
            print(f'Before filtering ligands: {len(results)}')
            results = filter_ligands(results, ligand_threshold)
            print(f'After filtering ligands with threshold of {ligand_threshold}: {len(results)}')
            print('Results len:', len(results), list(results.values())[0])
            print('Smiles len:', len(smiles))
            results_final, fails2 = pdb_ligands_vs_chembl_ligands(results, smiles, threshold=0.95, suffix=str(index))

            fails2 = pd.DataFrame(fails2, columns=['Target', 'Ligand_ID', 'Ligand_Smiles'])
            fails2.to_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, f'fail_smiles_pdb{index}.csv'))

            master_table = pd.read_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, file), index_col=0)

            zipped_results = zip_results_with_pdb(results_final, pdb_to_ligands, master_table)

            master_table = update_master_table(master_table, zipped_results)
            print('Starting primary PDB search')
            master_table = choose_primary_pdb_for_chembl(master_table)
            master_table.to_csv(os.path.join(COMPOUND_SAMPLING_FOLDER, f"targets_after_fingerprint_similarity{index}_tc{threshold}.csv"))
            evaluate(index, (threshold,))
            print(f'For index {index} - Targets having at least one ligand with TC > {threshold} '
                  f'to FP base of given target: {len(results_final)}')
    make_chembls_smiles_files(os.path.join(COMPOUND_SAMPLING_FOLDER,
                                           f'pdb_vs_chembl_filtered_ligands_tc{CHOSEN_TC_THRESHOLD}-{CHOSEN_LIGAND_LIMIT}.pkl'))
