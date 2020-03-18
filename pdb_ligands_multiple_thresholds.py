import rdkit
import pandas as pd
import pypdb
import pickle
from collections import Counter
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import rdkit.Chem as Chem
import os
import numpy as np
import time


def get_ligands(pdbs, suffix="", printing=False):
    results = dict()
    pdb_to_ligands = dict()
    s = 1
    for target in pdbs:
        if s % 100 == 0:
            # print('Currently checked ' + str(s) + ' targets.')
            pass
        ligands = []
        smiles = dict()
        ex = pdbs[target].split()
        # print(str(len(ex)) + ' PDB\'s found for target.')
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
                    # print('Failed to fetch ligands, retrying.')
                    if failes > 100:
                        # print(f'Can\'t download ligand {pdb_id}')
                        raise Exception('Connection to PDB server lost.')
                    continue
                break
            k += 1
            if printing:
                # print(k / len(ex), s / len(pdbs))
                pass
        results[target] = [ligands, smiles]
        s += 1
    pickle.dump(results, open(f"ligands{suffix}.pkl", "wb"))
    pickle.dump(pdb_to_ligands, open(f"pdb_to_ligands{suffix}.pkl", "wb"))
    return results, pdb_to_ligands


def count_ligands(results):
    res = []
    for target in results:
        ligs = results[target][0]
        res += ligs
    res = dict(Counter(res))
    return res


def make_fingerprints_base(targets, folder_name=False, index='', path='raw_data', prefix='CHEMBL25-single'):
    # print(f'Index for making fingerprints base {index}')
    files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and prefix in i]
    standard_csvs = [pd.read_csv(path + '/' + i, low_memory=False)[['Molecule', 'Target', 'Canonical Smiles']].dropna()
                     for i in files]
    kd = pd.read_csv('raw_data/CHEMBL25-kd.csv', delimiter=';', low_memory=False)
    kd.dropna()
    standard_csvs.append(kd)
    pdbs = list(targets['ChEMBL ID'])
    smile_dict = dict()
    fails = []
    for target in pdbs:
        smile_dict[target] = []
        p = set()
        s2 = dict()
        for standard in standard_csvs:
            target_standard_smiles = set(standard[standard['Target'] == target]['Canonical Smiles'])
            p.update(target_standard_smiles)
            s3 = standard[standard['Target'] == target][['Molecule', 'Canonical Smiles']]
            for index2, row in s3.iterrows():
                s2[row['Canonical Smiles']] = row['Molecule']
        p = list(p)
        k_inner = []
        for sm in p:
            try:
                chem_smiles = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(sm))
                k_inner.append((chem_smiles, sm, s2[sm]))
            except:
                # print('Failed to create Chem SMILES from: ', sm)
                fails.append((target, sm))
        smile_dict[target] = k_inner

    name_fp = f"fingerprints_chembl{index}.pkl"
    name_failes = f"fingerprints_chembl_fails{index}.pkl"

    if folder_name != False:
        name_fp = folder_name + '/' + name_fp
        name_failes = folder_name + '/' + name_failes
    # print(len(smile_dict))
    pickle.dump(smile_dict, open(name_fp, "wb"))
    pickle.dump(fails, open(name_failes, "wb"))
    # print(f'Saved fingerptints to {name_fp}')
    return smile_dict, fails


def pdb_ligands_vs_chembl_ligands(pdbs, fingerprints, threshold=0.95, to_folder=False, suffix=""):
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
                l2 = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(ligand_smile))
                for chembl_smile in target_chembl:
                    coef = DataStructs.FingerprintSimilarity(l2, chembl_smile[0])
                    if coef > threshold:
                        output[target][ligand_id].append((coef, chembl_smile[2]))
            except:
                # print(ligand_tuple)
                fails.append((target, ligand_tuple))
        for i in list(output[target].keys()):
            if len(output[target][i]) == 0:
                output[target].pop(i)
    for i in list(output.keys()):
        if len(output[i]) == 0:
            output.pop(i, None)
    fails = [(i[0], i[1][0], i[1][1]) for i in fails]

    if to_folder is not False:
        pickle.dump(output, open(f"{to_folder}/pdb_vs_chembl_filtered_ligands_tc{threshold}-{suffix}.pkl", "wb"))
        pickle.dump(fails, open(f"{to_folder}/fail_smiles_pdb_tc{threshold}{suffix}.pkl", "wb"))
    else:
        pickle.dump(output, open(f"pdb_vs_chembl_filtered_ligands_tc{threshold}-{suffix}.pkl", "wb"))
        pickle.dump(fails, open(f"fail_smiles_pdb_tc{threshold}{suffix}.pkl", "wb"))
    return output, fails


def filter_ligands(all_ligands, threshold=9):
    to_remove = set(list(dict((k, v) for k, v in count_ligands(all_ligands).items() if v > threshold).keys()))
    for target in list(all_ligands.keys()):
        keys = set(all_ligands[target][0])
        intersection = keys.intersection(to_remove)
        for lig_name in intersection:
            all_ligands[target][0].remove(lig_name)
            all_ligands[target][1].pop(lig_name, None)
        if len(all_ligands[target][0]) == 0:
            all_ligands.pop(target, None)
    return all_ligands


def evaluate(suffix, tcs_list=[0.95], folder='actives_number_sampling'):
    '''
    Make a matrix with comparison between runs on different parameters (how many proteins got through filtering).
    '''
    matrix = np.zeros((len(tcs_list), 1))
    path = folder + '/'
    for tc in range(len(tcs_list)):
        name = f'{path}pdb_vs_chembl_filtered_ligands_tc{tcs_list[tc]}-{suffix}.pkl'
        with open(name, 'rb') as handle:
            filtered_ligands = pickle.load(handle)
        good_targets = []
        for target in list(filtered_ligands.keys()):
            values = filtered_ligands[target]
            output = 0
            for val in values:
                if val[3] == 1:
                    output = 1
            if output == 0:
                good_targets.append(target)
        pickle.dump(good_targets, open(f"{path}targets_not_in_pdbbind_tc{tcs_list[tc]}-{suffix}.pkl", "wb"))
        with open(f"{path}targets_not_in_pdbbind_tc{tcs_list[tc]}-{suffix}.txt", 'w') as save_file:
            for i in good_targets:
                save_file.write(i + '\n')
        matrix[tc, 0] = len(good_targets)
    pickle.dump(matrix, open(f"{path}runs_comparison_matrix{suffix}.pkl", "wb"))
    pd.DataFrame(data=matrix, index=tcs_list).to_csv(f'{path}/runs_comparison_matrix{suffix}.csv')
    return matrix


def try_other_parameters(ligands, fingerprints, TCthreshold, folder=False, suffix=""):
    # print(f"Filtering {len(ligands)} ligands by parameters: Tanimoto Codefficient {TCthreshold}.")
    if folder is not False:
        if not os.path.exists(folder):
            os.makedirs(folder)

    results_filtered = ligands
    output, fails = pdb_ligands_vs_chembl_ligands(results_filtered, fingerprints, TCthreshold, to_folder=folder)
    # print(len(output))
    name = f"pdb_vs_chembl_filtered_ligands_tc{TCthreshold}.pkl"
    if folder is not False:
        name = f"{folder}/pdb_vs_chembl_filtered_ligands_tc{TCthreshold}.pkl"
    with open(name, 'rb') as handle:
        results_final = pickle.load(handle)

    with open('fail_smiles_chembl.pkl', 'rb') as handle:
        fails1 = pickle.load(handle)

    name = f"fail_smiles_pdb_tc{TCthreshold}.pkl"
    if folder is not False:
        name = f"{folder}/fail_smiles_pdb_tc{TCthreshold}.pkl"
    with open(name, 'rb') as handle:
        fails2 = pickle.load(handle)

    with open('pdb_to_ligands.pkl', 'rb') as handle:
        pdb_to_ligands = pickle.load(handle)

    fails1 = pd.DataFrame(fails1, columns=['Target', 'Ligand_Smiles'])  # 7 związków z chembl
    fails1.to_csv('fail_smiles_chembl.csv')

    name = f'fails_smiles_pdb_TC{TCthreshold}.csv'
    if folder is not False:
        name = f'{folder}/fails_smiles_pdb_TC{TCthreshold}.csv'
    fails2 = pd.DataFrame(fails2, columns=['Target', 'Ligand_ID', 'Ligand_Smiles'])
    fails2.to_csv(name)

    master_table = pd.read_csv('targets_without_coverage.csv', index_col=0)
    zipped_results = zip_results_with_pdb(results_final, pdb_to_ligands, master_table)
    pdbbind = read_pdbbind_ids('pdbbind_ids.txt')
    zipped_results = results_in_pdbbind(zipped_results, pdbbind)
    if folder is not False:
        pickle.dump(zipped_results, open(f"{folder}/pdb_vs_chembl_filtered_ligands_tc{TCthreshold}.pkl", "wb"))
    master_table = update_master_table(master_table, zipped_results, TCthreshold, folder=folder)
    if folder is not False:
        make_pdbbind_csv(zipped_results,
                         name=f'{folder}/targets_pdbs_ligands_tc{TCthreshold}.csv')
    else:
        make_pdbbind_csv(zipped_results, name=f'targets_pdbs_ligands_tc{TCthreshold}.csv')


def read_pdbbind_ids(file):
    ids = []
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '#':
                ids.append(line.split(' ')[0].upper())
    return ids


def zip_results_with_pdb(results_final, pdb_to_ligands, master_table):
    output = dict()
    for target in results_final:
        output[target] = []  # krotki (pdb_id, ligand, TC)
        try:
            possible_pdbs = list(master_table[master_table['ChEMBL ID'] == target]['PDB_entry'])[0].split()
        except:
            # print('No PDB entrys for: ', target)
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


def update_master_table(master_table, results, tc=0.95, folder=False, suffix=""):
    master_table = master_table[master_table['ChEMBL ID'].isin(list(results.keys()))].copy(deep=True)
    master_table['PDBs_number'] = 0
    master_table['Ligand_number'] = 0
    master_table['In_PDBBIND'] = 0
    for target in results:
        pdbs = set()
        ligs = set()
        in_pdbbind = set()
        for each in results[target]:
            pdbs.add(each[0])
            ligs.add(each[1])
            if each[3] == 1:
                in_pdbbind.add(each[0])
        master_table.at[target, 'PDBs_number'] = len(pdbs)
        master_table.at[target, 'Ligand_number'] = len(ligs)
        master_table.at[target, 'In_PDBBIND'] = len(in_pdbbind)
        new_pdbs = " ".join(list(pdbs))
        master_table.at[target, 'PDB_entry'] = new_pdbs
    if folder is not False:
        master_table.to_csv(f"{folder}/targets_after_fingerprint_similarity{suffix}_tc{tc}.csv")
    else:
        master_table.to_csv(f"targets_after_fingerprint_similarity{suffix}_tc{tc}.csv")
    return master_table


def results_in_pdbbind(zipped_results, pdbbind):
    output = dict()
    for target in zipped_results:
        output[target] = []
        for ligand in zipped_results[target]:
            target_pdb = ligand[0]
            target_ligand = ligand[1]
            target_tc = ligand[2]
            if target_pdb in pdbbind:
                output[target].append((target_pdb, target_ligand, target_tc, 1))
            else:
                output[target].append((target_pdb, target_ligand, target_tc, 0))
    return output


def make_pdbbind_csv(zipped_results, name='targets_pdbs_ligands_tc_pdbbind.csv'):
    table = pd.DataFrame(columns=['ChEMBL ID', 'PDB_ID', 'Ligand_ID', 'Tc_coef', 'Is_in_PDBBIND'])
    for target in zipped_results:
        for each in zipped_results[target]:
            table = table.append({'ChEMBL ID': target, 'PDB_ID': each[0], 'Ligand_ID': each[1],
                                  'Tc_coef': each[2], 'Is_in_PDBBIND': each[3]}, ignore_index=True)
    table.to_csv(name, index=False)


def sample_by_activity_conda(folder_with_masters, threshold=0.95):
    for index, file in enumerate([i for i in os.listdir(folder_with_masters) if 'targets_without_coverage' in i]):
        # print(index)
        data = pd.read_csv(f'{folder_with_masters}/{file}', index_col=0)
        pdbs = dict(data['PDB_entry'])
        try:
            with open(f'ligands{index}.pkl', 'rb') as handle:
                results = pickle.load(handle)
            with open(f'pdb_to_ligands{index}.pkl', 'rb') as handle:
                pdb_to_ligands = pickle.load(handle)
        except:
            results, pdb_to_ligands = get_ligands(pdbs, suffix=str(index), printing=False)  # download ligands
            print(f'For index {index} - Number of targets with ligands: {len(results)} , previously 808')

        try:
            with open(f'{folder_with_masters}/fingerprints_chembl.pkl', 'rb') as handle:
                smiles = pickle.load(handle)
                print('Fingerprints found and loaded!')
        except:
            print('Preparing fingerprint base...')
            smiles, fails1 = make_fingerprints_base(targets=pd.read_csv(f'master_table.csv', index_col=0),
                                                    folder_name=folder_with_masters)

        print('Results len:', len(results), list(results.values())[0])
        print('Smiles len:', len(smiles))
        results_final, fails2 = pdb_ligands_vs_chembl_ligands(results, smiles, threshold=0.95,
                                                              to_folder=folder_with_masters, suffix=str(index))
        print(f'For index {index} - Targety mające przynajmniej jeden ligand o punktacji TC >0.95 do bazy FP dla danego '
              f'targetu: {len(results_final)} , previously 390')

        fails2 = pd.DataFrame(fails2, columns=['Target', 'Ligand_ID', 'Ligand_Smiles'])
        fails2.to_csv(f'{folder_with_masters}/fail_smiles_pdb{index}.csv')

        master_table = pd.read_csv(f'{folder_with_masters}/{file}', index_col=0)

        zipped_results = zip_results_with_pdb(results_final, pdb_to_ligands, master_table)
        pdbbind = read_pdbbind_ids('pdbbind_ids.txt')
        zipped_results = results_in_pdbbind(zipped_results, pdbbind)
        pickle.dump(zipped_results,
                    open(f"{folder_with_masters}/pdb_vs_chembl_filtered_ligands_tc{threshold}-{index}.pkl", "wb"))
        master_table = update_master_table(master_table, zipped_results, folder=folder_with_masters, suffix=str(index))
        make_pdbbind_csv(zipped_results, name=f'{folder_with_masters}/targets_pdbs_ligands_tc_0.95_pdbbind{index}.csv')

        master_table = master_table[master_table['In_PDBBIND'] == 0]
        master_table.to_csv(f'master_table_not_in_pdbbind_{index}.csv')
        print(f'For index {index} - Number of targets that are not in PDBBIND: {len(master_table)}')

        evaluate(index, [threshold], folder=folder_with_masters)

sample_by_activity_conda('actives_number_sampling')