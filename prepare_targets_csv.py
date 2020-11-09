import os
import numpy as np
from paths_and_settings import *
from utils import create_folder_if_not_existent
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class Standard:
    def __init__(self):
        self.std = ''
        self.relation = ''
        self.value = 0
        self.parent = ''
        self.activity = None


def filter_csv_by_ids(filtered, loaded_csv, target='ChEMBL ID'):
    return loaded_csv[loaded_csv[target].isin(filtered)]


def filter_csv_by_units(desired_units, loaded_csv, target='Standard Units'):
    return loaded_csv[loaded_csv[target].isin(desired_units)]


def add_compounds_with_std(main_csv, loaded_csv, standard_type='Ki', std_id='Target'):
    main_csv['{}_all'.format(standard_type)] = 0
    for index, row in loaded_csv.iterrows():
        target_id = row[std_id]
        main_csv.at[target_id, '{}_all'.format(standard_type)] += 1
    return main_csv


def mol_weight_from_smiles(smile):
    x = Chem.MolFromSmiles(smile)
    return rdMolDescriptors.CalcExactMolWt(x)  # return grams per mol


def unit_converter(value, value_type, series):
    try:
        value = value.replace(',', '.')
    except:
        pass
    if value == np.nan:
        return np.nan
    elif value_type == 'nM' or value_type == '/nM':
        value = float(value) / 1000
    elif value_type == 'ug.mL-1' or value_type == 'ug ml-1':
        mol_weight = mol_weight_from_smiles(series['Canonical Smiles'])
        value = (float(value) * 1000) / mol_weight
    return float(value)


def check_activity_value_relation(main_csv, loaded_csv, standard_type, std_id='Target',
                                  units='Standard Units', value='Standard Value', relation='Standard Relation'):
    main_csv['{}_active'.format(standard_type)] = 0
    main_csv['{}_inactive'.format(standard_type)] = 0
    for index, row in loaded_csv.iterrows():
        target_id = row[std_id]
        rel = row[relation]  # =, >, <
        std_value = row[value]
        try:
            std_value = std_value.replace(',', '.')
        except AttributeError:
            pass
        if std_value == np.nan:
            main_csv.at[target_id, '{}_inactive'.format(standard_type)] += 1
            continue
        elif row[units] == 'nM' or row[units] == '/nM':
            std_value = float(std_value) / 1000
        elif row[units] == 'ug.mL-1' or row[units] == 'ug ml-1':
            mol_weight = mol_weight_from_smiles(row['Canonical Smiles'])
            std_value = (float(std_value) * 1000) / mol_weight

        std_value = float(std_value)
        if std_value < ACTIVITY_VALUE_THRESHOLD and rel == '=':
            main_csv.at[target_id, '{}_active'.format(standard_type)] += 1
            continue
        elif std_value >= 10 and rel == '=':
            main_csv.at[target_id, '{}_inactive'.format(standard_type)] += 1
            continue
    return main_csv


def make_standard_csv(loaded_csv, standard_type, compounds='Molecule', std_id='Target', units='Standard Units',
                      value='Standard Value', relation='Standard Relation', smile='Canonical Smiles'):
    """
    Make a CSV with rows ID_compound, ID_target_CHEMBL, SMILES, ACTIVE (value, type)/INACTIVE (value, type)
    """

    header = ['ID_compound', 'ID_target_CHEMBL', 'Target_SMILES', 'Activity', 'Value', 'Relation', 'Standard']
    create_folder_if_not_existent(STANDARDS_FOLDER)
    with open(os.path.join(STANDARDS_FOLDER, f'compounds_{standard_type}.csv'), 'w') as file:
        file.write(','.join(header))
        file.write('\n')
        for _, row in loaded_csv.iterrows():
            compound_id = row[compounds]
            target_id = row[std_id]
            target_smile = row[smile]
            rel = row[relation]  # =, >, <
            std_value = row[value]
            to_save = [compound_id, target_id, target_smile]
            if type(std_value) is str:
                std_value = std_value.replace(',', '.')
            if std_value is np.nan:
                to_save.append('Inactive')
                to_save.append('')
                to_save.append('')
                to_save.append(str(standard_type))
                for ix, i in enumerate(to_save):
                    if i == np.nan:
                        to_save[ix] = ''
                    else:
                        to_save[ix] = str(to_save[ix])
                file.write(','.join(to_save))
                file.write('\n')
                continue
            elif row[units] == 'nM' or row[units] == '/nM':
                std_value = float(std_value) / 1000
            elif row[units] == 'ug.mL-1' or row[units] == 'ug ml-1':
                mol_weight = mol_weight_from_smiles(row['Canonical Smiles'])
                std_value = float(std_value) * 1000 / mol_weight

            std_value = float(std_value)
            if std_value < ACTIVITY_VALUE_THRESHOLD and rel == '=':
                to_save.append('Active')
            elif std_value < ACTIVITY_VALUE_THRESHOLD and rel != '=':
                continue
            else:
                to_save.append('Inctive')

            to_save.append(str(std_value))
            to_save.append(str(rel))
            to_save.append(str(standard_type))
            for ix, i in enumerate(to_save):
                if i == np.nan or i == 'nan':
                    to_save[ix] = ''
                else:
                    to_save[ix] = str(to_save[ix])
            file.write(','.join(to_save))
            file.write('\n')


def sum_act_inact(main_csv, standards=('Ki', 'Kd', 'IC50')):
    for _, row in main_csv.iterrows():
        target_id = row['ChEMBL ID']
        summed_active = 0
        summed_inactive = 0
        for i in standards:
            summed_active += row['{}_active'.format(i)]
            summed_inactive += row['{}_inactive'.format(i)]
        main_csv.at[target_id, 'Active_compounds'] = summed_active
        main_csv.at[target_id, 'Inactive_compounds'] = summed_inactive
    return main_csv


def decois_uniID_from_folder(dir_name, extension, substring1='Name', substring2='Uniprot_ID',
                             standards=('IC50', 'Ki', 'Kd')):
    decoy_dir = os.path.join(dir_name, 'decoys')
    ligands_dir = os.path.join(dir_name, 'ligands')
    name_dict = dict()
    uniprot_ids = set()
    uni_dict = dict()
    uni_dict2 = dict()
    for item in os.listdir(ligands_dir):  # loop through items in dir
        if item.endswith(extension):  # check for ".gz" extension
            name = item.split(extension)[0]
            file_name = os.path.join(ligands_dir, item)  # get full path of files
            with open(file_name, 'r') as sdf_file:
                s = 0
                for line in sdf_file:
                    if s == 1:
                        name_temp = line.strip('\n')
                        try:
                            uni_dict[name_temp]
                        except KeyError:
                            uni_dict[name_temp] = []
                        uni_dict2[name_temp] = []
                        s = 0
                    if s == 2:
                        uni_temp = line.strip('\n').split(' ')
                        for i in uni_temp:
                            uniprot_ids.add(i)
                            try:
                                uni_dict2[name_temp].append(i)
                            except KeyError:
                                uni_dict2[name_temp] = i
                            if new_std.parent == '':
                                new_std.parent = i
                            else:
                                new_std.parent = ' '.join([new_std.parent, i])
                            lll = i
                        s = 0
                    if s == 3:
                        try:
                            new_std.value = float(line.strip('\n')) / 1000
                            new_std.relation = '='
                            uni_dict[name_temp].append(new_std)
                        except:
                            if '>' in line:
                                new_std.value = (float(line.strip('\n').strip('>')) / 1000) + 1
                                if new_std.value < ACTIVITY_VALUE_THRESHOLD:
                                    new_std.activity = 1
                                else:
                                    new_std.activity = 0
                                new_std.relation = '>'
                                uni_dict[name_temp].append(new_std)
                        s = 0
                    if substring2 in line:
                        s = 2
                    if substring1 in line:
                        s = 1
                    for std in standards:
                        if std in line:
                            new_std = Standard()
                            s = 3
                            new_std.std = std
            decoys_number = 0
            for pf in os.listdir(decoy_dir):  # loop through items in dir
                if name in pf:
                    decoys_number = 0
                    with open(os.path.join(decoy_dir, pf), 'r') as sdf_file:
                        for line in sdf_file:
                            if 'Name' in line:
                                decoys_number += 1
                    break
            name_dict[name] = [lll, decoys_number]
    return uni_dict2, name_dict, uni_dict


def decois_active_and_threshold(main_csv, decois_ligands, decois_ligands_standard):
    main_csv['DEKOIS_actives'] = 0
    main_csv['DEKOIS_actives_threshold'] = 0
    main_uni_ids = list(main_csv['UniProt Accession'])
    for key in decois_ligands:
        if decois_ligands[key] in main_uni_ids:
            for index, row in main_csv.iterrows():
                if row['UniProt Accession'] == decois_ligands[key]:
                    main_csv.at[index, 'DEKOIS_actives'] = main_csv.at[index, 'DEKOIS_actives'] + 1
                    select_standards = decois_ligands_standard[key]
                    for std in select_standards:
                        if std.value < ACTIVITY_VALUE_THRESHOLD:
                            main_csv.at[index, 'DEKOIS_actives_threshold'] = main_csv.at[
                                                                                 index, 'DEKOIS_actives_threshold'] + 1
                            break
                    break
    return main_csv


def decois_decoys_number_to_master(main_csv, decois_name_dict, decois_ligands, decois_ligands_standard):
    main_csv['DEKOIS_ID'] = ''
    main_csv = decois_active_and_threshold(main_csv, decois_ligands, decois_ligands_standard)
    main_csv['DEKOIS_decoys'] = 0
    main_uni_ids = list(main_csv['UniProt Accession'])
    for key in decois_name_dict:
        if decois_name_dict[key][0] in main_uni_ids:
            for index, row in main_csv.iterrows():
                if row['UniProt Accession'] == decois_name_dict[key][0]:
                    main_csv.at[index, 'DEKOIS_ID'] = key
                    main_csv.at[index, 'DEKOIS_decoys'] = decois_name_dict[key][1]
    return main_csv


def check_if_dude_in_master(main_csv, dude_path=DUDE_PATH,
                            standards=('Kd =', 'Ki =', 'IC50 =', 'Kd <', 'Ki <', 'IC50 <', 'Kd >', 'Ki >', 'IC50 >')):
    k = ['actives_nM_chembl.ism', 'inactives_nM_chembl.ism']
    k2 = ['Active_in_DUDE', 'Active_in_DUDE_threshold', 'Inactive_in_DUDE', 'Inactive_in_DUDE_threshold']
    bfolders = [f.path for f in os.scandir(dude_path) if f.is_dir()]
    main_csv['DUDE_ID'] = ''
    for i in k2:
        main_csv[i] = 0
    for dude in bfolders:
        dude_name = dude.split('/')[-1]
        with open(os.path.join(dude_path, dude_name, 'uniprot.txt')) as f:
            lineList = [l.strip('\n') for l in f.readlines()]
        for index, row in main_csv.iterrows():
            uni_id = row['UniProt Accession']
            if uni_id in lineList:
                if len(main_csv.at[index, 'DUDE_ID']) > 0:
                    main_csv.at[index, 'DUDE_ID'] = main_csv.at[index, 'DUDE_ID'] + ' {}'.format(dude_name)
                else:
                    main_csv.at[index, 'DUDE_ID'] = '{}'.format(dude_name)
                for file in range(len(k)):
                    with open(os.path.join(dude_path, dude_name, k[file])) as f:
                        lineList = [l.strip('\n') for l in f.readlines()]
                        if file == 0:
                            main_csv.at[index, k2[0]] = main_csv.at[index, k2[0]] + len(lineList)
                            for line in lineList:
                                for std in standards:
                                    if std in line:
                                        value = float(line.split(std)[1].split(' ')[1]) / 1000
                                        if value < ACTIVITY_VALUE_THRESHOLD:
                                            main_csv.at[index, k2[1]] = main_csv.at[index, k2[1]] + 1
                                        break
                        if file == 1:
                            main_csv.at[index, k2[2]] = main_csv.at[index, k2[2]] + len(lineList)
                            for line in lineList:
                                for std in standards:
                                    if std in line:
                                        value = float(line.split(std)[1].split(' ')[1]) / 1000
                                        if value >= 10:
                                            main_csv.at[index, k2[3]] = main_csv.at[index, k2[3]] + 1
                                        break
    return main_csv


# 'Document Year'  'Standard Value', 'Standard Units', 'Standard Relation'
def duplicates_parser(read_csv):
    duplicates = read_csv[read_csv.duplicated(keep=False, subset=['Target', 'Molecule'])]
    ids = set(duplicates['Target'])
    to_drop = []
    for i in ids:
        dup_a = duplicates[duplicates['Target'] == i]
        dup_mol_id = set(dup_a['Molecule'])
        for j in dup_mol_id:
            dup = dup_a[dup_a['Molecule'] == j]
            dup_values = dup[['Molecule', 'Target']]
            dup_ids = list(dup.index.tolist())
            if len(dup_values) > 1:
                main_index = dup_ids[0]
                main_values = {i: dup.loc[main_index][i] for i in dup.loc[main_index].keys()}
                for ix in range(1, len(dup_ids)):
                    diff = (unit_converter(main_values['Standard Value'], main_values['Standard Units'],
                                           dup.loc[main_index]) + 1e-10) / (
                                       unit_converter(dup.loc[dup_ids[ix]]['Standard Value'],
                                                      dup.loc[dup_ids[ix]]['Standard Units'],
                                                      dup.loc[dup_ids[ix]]) + 1e-10)
                    if diff != np.nan or diff < 0.8 or diff > 1.2:
                        main_year = main_values['Document Year']
                        new_year = dup.loc[dup_ids[ix]]['Document Year']
                        if new_year > main_year:
                            main_index = dup_ids[ix]
                            main_values = {i: dup.loc[main_index][i] for i in dup.loc[main_index].keys()}
                [to_drop.append(i) for i in dup_ids]
    read_csv = read_csv.drop(to_drop)
    return read_csv
