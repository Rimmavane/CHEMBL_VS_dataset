import pandas as pd
import pickle
import os


def make_chembls_smiles_files(path_to_fitered_targets, smiles_folder='standards_csv', prefix='compounds', folder='chembl_smiles'):
    """
    Make SMILES file for each chembl id from selected file in pickle format.
    :param path_to_fitered_targets: Pickled list of ids.
    :param smiles_folder: Where to find smiles csv's.
    :param prefix: Where to find csv's with smiles.
    :param folder: Where to save smiles files.
    :return:
    """
    if not os.path.exists('chembl_smiles'):
        os.makedirs('chembl_smiles')
    with open(path_to_fitered_targets, 'rb') as handle:
        filtered_targets = pickle.load(handle)
    files = [i for i in os.listdir(smiles_folder) if os.path.isfile(os.path.join(smiles_folder, i)) and prefix in i] # ZAŁADOWANE TARGET CHEMBL ID, ZNALEŹĆ ICH SMILES Z CHEMBL (?) CZY LIGANDOWE (??)
    standard_csvs = [pd.read_csv(smiles_folder + '/' + i, low_memory=False).dropna() for i in files]
    smile_dict = dict()
    for chembl_id in filtered_targets:
        smile_dict[chembl_id] = dict()
        with open(f'{folder}/{chembl_id}_all.smi', 'w') as save_file:
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_TARGET_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        smile = record[1]['Target_SMILES']
                        chembl_ligand = record[1]['ID_compound']
                        save_file.write(smile + '\t' + chembl_ligand + '\n')
                        smile_dict[chembl_id][chembl_ligand] = smile
                except KeyError:
                    pass
        with open(f'{folder}/{chembl_id}_active.smi', 'w') as save_file:
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_TARGET_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        if record[1]['Activity'] == 'Active':
                            smile = record[1]['Target_SMILES']
                            chembl_ligand = record[1]['ID_compound']
                            save_file.write(smile + '\t' + chembl_ligand + '\n')
                except KeyError:
                    pass
        with open(f'{folder}/{chembl_id}_inactive.smi', 'w') as save_file:
            for standard in standard_csvs:
                try:
                    good_records = standard[standard['ID_TARGET_CHEMBL'] == chembl_id]  # ['ID_compound', 'Target_SMILES']
                    for record in good_records.iterrows():
                        if record[1]['Activity'] != 'Active':
                            smile = record[1]['Target_SMILES']
                            chembl_ligand = record[1]['ID_compound']
                            save_file.write(smile + '\t' + chembl_ligand + '\n')
                except KeyError:
                    pass
    return smile_dict


smiles = make_chembls_smiles_files('parameters_try/targets_not_in_pdbbind_tc0.95.pkl')
