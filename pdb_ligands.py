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


def get_ligands(pdbs):
    results = dict()
    pdb_to_ligands = dict()
    s = 1
    for target in pdbs:
        if s % 100 == 0:
            print('Currently checked ' + str(s) + ' targets.')
        ligands = []
        smiles = dict()
        ex = pdbs[target].split()
        print(str(len(ex)) + ' PDB\'s found for target.')
        k = 0
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
                except:
                    print('Failed to fetch ligands, retrying.')
                    continue
                break
            k += 1
            print(k / len(ex), s / len(pdbs))
        results[target] = [ligands, smiles]
        s += 1
    pickle.dump(results, open("ligands.pkl", "wb"))
    pickle.dump(results, open("pdb_to_ligands.pkl", "wb"))
    return results, pdb_to_ligands


def count_ligands(results):
    res = []
    for target in results:
        ligs = results[target][0]
        res += ligs
    res = dict(Counter(res))
    return res


def make_fingerprints_base(targets, path='raw_data', prefix='CHEMBL25-single'):
    files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and prefix in i]
    standard_csvs = [pd.read_csv(path + '/' + i, low_memory=False)[['Molecule', 'Target', 'Canonical Smiles']].dropna() for i in files]
    kd = pd.read_csv('raw_data/CHEMBL25-kd.csv', delimiter=';', low_memory=False)
    kd['Target'] = kd['Target ChEMBL ID']
    kd['Molecule'] = kd['Molecule ChEMBL ID']
    kd['Canonical Smiles'] = kd['Smiles']
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
            for index, row in s3.iterrows():
                s2[row['Canonical Smiles']] = row['Molecule']
        p = list(p)
        k_inner = []
        for sm in p:
            try:
                chem_smiles = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(sm))
                k_inner.append((chem_smiles, sm, s2[sm]))
            except:
                print(sm)
                fails.append((target, sm))
        smile_dict[target] = k_inner
    pickle.dump(smile_dict, open("fingerprints_chembl.pkl", "wb"))
    pickle.dump(fails, open("fail_smiles_chembl.pkl", "wb"))
    return smile_dict, fails


def pdb_ligands_vs_chembl_ligands(pdbs, fingerprints, threshold=0.95, to_folder=False):
    targets = list(pdbs.keys())
    output = dict()
    fails = []
    for target in targets:
        output[target] = dict()
        target_smiles = list(pdbs[target][1].items())  #smiles ligandów z PDB
        target_chembl = fingerprints[target]
        for ligand_tuple in target_smiles:
            try:
                ligand_id = ligand_tuple[0] #id ligandu
                ligand_smile = ligand_tuple[1]  #smiles ligandu
                output[target][ligand_id] = []
                l2 = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(ligand_smile))
                for chembl_smile in target_chembl:
                    coef = DataStructs.FingerprintSimilarity(l2, chembl_smile[0])
                    if coef > threshold:
                        output[target][ligand_id].append((coef, chembl_smile[2]))
            except:
                print(ligand_tuple)
                fails.append((target, ligand_tuple))
        for i in list(output[target].keys()):
            if len(output[target][i]) == 0:
                output[target].pop(i)
    for i in list(output.keys()):
        if len(output[i]) == 0:
            output.pop(i, None)
    fails = [(i[0], i[1][0], i[1][1]) for i in fails]

    if to_folder is not False:
        pickle.dump(output, open(f"{to_folder}/pdb_vs_chembl_filtered_ligands_tc{threshold}.pkl", "wb"))
        pickle.dump(fails, open(f"{to_folder}/fail_smiles_pdb_tc{threshold}.pkl", "wb"))
    else:
        pickle.dump(output, open(f"pdb_vs_chembl_filtered_ligands_tc{threshold}.pkl", "wb"))
        pickle.dump(fails, open(f"fail_smiles_pdb_tc{threshold}.pkl", "wb"))
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


def evaluate(tcs_list=[0.95, 0.9, 0.85, 0.75, 0.7, 0.65, 0.6], folder='parameters_try'):
    '''
    Make a matrix with comparison between runs on different parameters (how many proteins got through filtering).
    '''
    matrix = np.zeros((len(tcs_list), 1))
    path = folder+'/'
    for tc in range(len(tcs_list)):
        name = f'pdb_vs_chembl_filtered_ligands_tc{tcs_list[tc]}.pkl'
        with open(path+name, 'rb') as handle:
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
        pickle.dump(good_targets, open(f"{path}targets_not_in_pdbbind_tc{tcs_list[tc]}.pkl", "wb"))
        with open(f"{path}targets_not_in_pdbbind_tc{tcs_list[tc]}.txt", 'w') as save_file:
            for i in good_targets:
                save_file.write(i+'\n')
        matrix[tc, 0] = len(good_targets)
    pickle.dump(matrix, open(f"runs_comparison_matrix.pkl", "wb"))
    pd.DataFrame(data=matrix, index=tcs_list).to_csv('comparison_between_runs.csv')
    return matrix


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
            print(best_pdb)
            main_table.at[index, 'main_PDB_structure'] = best_pdb
    main_table.to_csv(csv_path)


def try_other_parameters(ligands, fingerprints, TCthreshold, folder=False):
    print(f"Filtering {len(ligands)} ligands by parameters: Tanimoto Codefficient {TCthreshold}.")
    if folder is not False:
        if not os.path.exists(folder):
            os.makedirs(folder)

    results_filtered = ligands
    output, fails = pdb_ligands_vs_chembl_ligands(results_filtered, fingerprints, TCthreshold, to_folder=folder)
    print(len(output))
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
    choose_primary_pdb_for_chembl('targets_without_coverage.csv')

    master_table = pd.read_csv('targets_without_coverage.csv', index_col=0)
    zipped_results = zip_results_with_pdb(results_final, pdb_to_ligands, master_table)
    pdbbind = read_pdbbind_ids('pdbbind_ids.txt')
    zipped_results = results_in_pdbbind(zipped_results, pdbbind)
    if folder is not False:
        pickle.dump(zipped_results, open(f"{folder}/pdb_vs_chembl_filtered_ligands_tc{TCthreshold}.pkl", "wb"))
    master_table = update_master_table(master_table, zipped_results, TCthreshold, folder=folder)
    if folder is not False:
        make_pdbbind_csv(master_table,
                         name=f'{folder}/targets_pdbs_ligands_tc{TCthreshold}.csv')
    else:
        make_pdbbind_csv(master_table, name=f'targets_pdbs_ligands_tc{TCthreshold}.csv')


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
        output[target] = []     #krotki (pdb_id, ligand, TC)
        possible_pdbs = list(master_table[master_table['ChEMBL ID'] == target]['PDB_entry'])[0].split()
        target_ligands = list(results_final[target].keys())     #ligand
        for ligand in target_ligands:
            for pdb_id in possible_pdbs:
                try:
                    if ligand in pdb_to_ligands[pdb_id]:
                        output[target].append((pdb_id, ligand, results_final[target][ligand]))
                except:
                    pass
    return output


def update_master_table(master_table, results, tc=0.95, folder=False):
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
        master_table.to_csv(f"{folder}/targets_after_fingerprint_similarity_tc{tc}.csv")
    else:
        master_table.to_csv(f"targets_after_fingerprint_similarity_tc{tc}.csv")
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


k = pd.read_csv('targets_without_coverage.csv', index_col=0)
pdbs = dict(k['PDB_entry'])
#results, pdb_to_ligands = get_ligands(pdbs) #download ligands

with open('ligands.pkl', 'rb') as handle:
    results = pickle.load(handle)       #808 targetów

smiles, fails1 = make_fingerprints_base(k)

#with open('fingerprints_chembl.pkl', 'rb') as handle:
#    smiles = pickle.load(handle)

print('results', len(results))
#results_filtered = filter_ligands(results, 9)         #626 targetów - odfiltrowano ligandy występujące >9 razy i targety bez ligandów

results_final, fails2 = pdb_ligands_vs_chembl_ligands(results, smiles)  #390 targetów - tylko targety mające przynajmniej jeden ligand o punktacji TC >0.95 do bazy FP dla danego targetu
print('results_final', len(results_final))

#with open('pdb_vs_chembl_filtered_ligands.pkl', 'rb') as handle:
#    results_final = pickle.load(handle)

#with open('fail_smiles_chembl.pkl', 'rb') as handle:
#    fails1 = pickle.load(handle)

#with open('fail_smiles_pdb.pkl', 'rb') as handle:
#    fails2 = pickle.load(handle)

with open('pdb_to_ligands.pkl', 'rb') as handle:
    pdb_to_ligands = pickle.load(handle)

fails1 = pd.DataFrame(fails1, columns=['Target', 'Ligand_Smiles'])  #7 związków z chembl
fails1.to_csv('fail_smiles_chembl.csv')

fails2 = pd.DataFrame(fails2, columns=['Target', 'Ligand_ID',  'Ligand_Smiles'])
fails2.to_csv('fail_smiles_pdb.csv')

master_table = pd.read_csv('targets_without_coverage.csv', index_col=0)
zipped_results = zip_results_with_pdb(results_final, pdb_to_ligands, master_table)
pdbbind = read_pdbbind_ids('pdbbind_ids.txt')

zipped_results = results_in_pdbbind(zipped_results, pdbbind)
master_table = update_master_table(master_table, zipped_results)
master_table.to_csv('master_table_updated.csv')
make_pdbbind_csv(zipped_results)
'''
targets_without_coverage
target - zipped_results.keys
description - master_table['ChEMBL ID'] --> 'Name'
liczba aktywnych - master_table['ChEMBL ID'] --> 'Active_compounds'
liczba nieaktywnych - master_table['ChEMBL ID'] --> 'Inactive_compounds'
lista chembl_ligand
lista PDB target-active ligand, 
jakie PDB jest w pdbbind - wynik results_in_pdbbind

'''

tcs = [0.95, 0.9, 0.85, 0.75, 0.7, 0.65, 0.6]
for y in tcs:
    with open('ligands.pkl', 'rb') as handle:
        results = pickle.load(handle)  # 808 targetów
        print('results pickled', len(results))
    with open('fingerprints_chembl.pkl', 'rb') as handle:
        smiles = pickle.load(handle)
    try_other_parameters(ligands=results, fingerprints=smiles, TCthreshold=y, folder='parameters_try')
evaluate(tcs)
