import pandas as pd
import os
from utils import log, load_smiles
from paths_and_settings import *


def load_activities():
    activities = []
    activities_paths = [join(STANDARDS_FOLDER, f) for f in os.listdir(STANDARDS_FOLDER)
                        if (os.path.isfile(join(STANDARDS_FOLDER, f)) and 'compound' in f)]
    for p in activities_paths:
        data = pd.read_csv(p)
        data['Value'] = pd.to_numeric(data['Value'])
        data = data[data['Relation'] == '=']
        data = data[data['Value'] < ACTIVITY_VALUE_THRESHOLD]
        activities.append(data)
    output = pd.concat(activities, ignore_index=True)
    return output


def get_ligands_best_activity(ligand_id: str, acts_dataframe: pd.DataFrame):
    best = min(acts_dataframe[acts_dataframe['ID_compound'] == ligand_id]['Value'].to_list())
    return best


def choose_best_actives(target_id: str, path_to_analog_matrix: str, acts_dataframe: pd.DataFrame, tc_threshold=ACTIVES_TC_SIMILARITY_THRESHOLD):
    chembl_actives_matrix = pd.read_csv(path_to_analog_matrix, index_col=0).apply(pd.to_numeric).to_dict(orient='index')
    log(f'Before filtering best actives for {target_id}: {len(chembl_actives_matrix)}')
    target_act_dataframe = acts_dataframe[acts_dataframe['ID_target_CHEMBL'] == target_id]
    for i in range(len(chembl_actives_matrix)):
        bail = False
        to_drop = None
        for query_id_1 in chembl_actives_matrix:
            for subject_id_2 in chembl_actives_matrix[query_id_1]:
                if query_id_1 != subject_id_2:
                    tc = max(chembl_actives_matrix[query_id_1][subject_id_2], chembl_actives_matrix[subject_id_2][query_id_1])
                    if tc >= tc_threshold:
                        best1 = get_ligands_best_activity(query_id_1, target_act_dataframe)
                        best2 = get_ligands_best_activity(subject_id_2, target_act_dataframe)
                        if best1 >= best2:
                            to_drop = subject_id_2
                        else:
                            to_drop = query_id_1
                        bail = True
                    if bail:
                        break
            if bail:
                break
        if to_drop is not None:
            del chembl_actives_matrix[to_drop]
            for k in chembl_actives_matrix:
                del chembl_actives_matrix[k][to_drop]
        else:
            break
    log(f'After filtering best actives for {target_id}: {len(chembl_actives_matrix)}')
    return list(chembl_actives_matrix.keys())


def filter_actives_smiles_file(target_id: str, best_actives: list):
    smiles = load_smiles(os.path.join(CHEMBL_SMILES_FOLDER, f'{target_id}_active.smi'))
    to_del = []
    for i in smiles:
        if i not in best_actives:
            to_del.append(i)
    for i in to_del:
        del smiles[i]
    with open(os.path.join(CHEMBL_SMILES_FOLDER, f'{target_id}_filtered_active.smi'), 'w') as handle:
        for tup in smiles.items():
            handle.write(f'{tup[1]}\t{tup[0]}\n')
