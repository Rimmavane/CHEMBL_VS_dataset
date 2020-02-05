import os
import pandas as pd
from time import time
from itertools import islice
import argparse
import pickle
import distutils.util

from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
import rdkit.Chem.Scaffolds.MurckoScaffold as murco
from rdkit.Chem.Fingerprints import FingerprintMols
from joblib import Parallel, delayed


def log(file_path='', string=''):
    if file_path != '':
        os.popen(f"echo {string} >> {file_path} &")
    else:
        print(string)


def get_files_with_suffix(folder_path, suffix='_active.smi'):
    """
    Find all files with given suffix in specified folder, without recursion
    :param folder_path: Folder we would like to search in
    :param suffix: Suffix of files we want to find
    :return: List of paths to files
    """
    files = [folder_path + '/' + i for i in os.listdir(folder_path) if
             os.path.isfile(os.path.join(folder_path, i)) and suffix in i]
    return files


class Decoy_finder:
    def __init__(self):
        pass

    def add_descriptors(self, targets_dict):
        for target in targets_dict:
            desc = self.make_descriptors_for_smile(targets_dict[target]['smile'])
            targets_dict[target].update(desc)
        return targets_dict

    def load_smiles_in_batches(self, file_path, skip_lines=0, batch_size=10000):
        with open(file_path, 'r') as handle:
            for _ in range(skip_lines):
                next(handle)
            i = 0
            smiles = dict()
            for line in handle:
                smile, smile_id = line.strip().split()
                smile_tuple = (smile, smile_id)
                smiles[smile_tuple[1]] = {'smile': smile_tuple[0], 'fingerprint': self.make_fingerprint(smile_tuple[0])}
                i += 1
                if i == batch_size:
                    yield smiles
                    i = 0
                    smiles = dict()
                    continue
            yield smiles

    def make_fingerprint(self, smile):
        """
        Make MACCS fingerprint for given smile
        :param smile: Compound smile
        :return: MACCS fingerprint
        """
        try:
            x = Chem.MolFromSmiles(smile)
            x = MACCSkeys.GenMACCSKeys(x)
            return x
        except:
            print(f'COULDNT MAKE FINGERPRINT FROM {smile}')
            return None

    def make_descriptors_for_smile(self, smile):
        """
        Make physical descriptors:
        - molecular weight,
        - the number of rotational bonds,
        - total hydrogen bond donors (HBDs),
        - total hydrogen bond acceptors (HBAs),
        - octanol–water partition coefficient (log P)
        :param smile: SMILE representation of compound
        :return: Dictionary with physical descriptors of smile
        """
        descriptors = dict()
        x = Chem.MolFromSmiles(smile)
        descriptors['HBD'] = rdMolDescriptors.CalcNumHBD(x)
        descriptors['HBA'] = rdMolDescriptors.CalcNumLipinskiHBA(x)
        descriptors['rotates'] = rdMolDescriptors.CalcNumRotatableBonds(x)
        descriptors['weight'] = rdMolDescriptors.CalcExactMolWt(x)
        descriptors['logp'] = rdMolDescriptors.CalcCrippenDescriptors(x)[0]
        descriptors['murco_scaffold'] = FingerprintMols.FingerprintMol(murco.GetScaffoldForMol(x))
        return descriptors


class Chembl_loader(Decoy_finder):
    def __init__(self, chembl_csv, chembl_smi_path, decoy_folder_path, th_dicts,
                 max_decoys=50, log_file='', ignore_pickled=True):
        super().__init__()
        self.log_file = log_file
        self.n_jobs = os.cpu_count() - 1
        self.global_start_time = time()
        self.chembl_smiles_path = chembl_smi_path
        self.decoy_base_folder_path = decoy_folder_path
        self.threshold_dicts = th_dicts
        self.max_decoys_for_active = max_decoys
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'): pass
        log(self.log_file, f"Availible cores: {self.n_jobs}")

        log(self.log_file, 'Loading chembl ids...')
        if type(chembl_csv) != list:
            self.targets_ids = self.get_chembl_ids_from_csv(chembl_csv)  # load chembl ids
        else:
            self.targets_ids = chembl_csv

        ### TEST SUBJECT ##
        self.targets_ids = ['CHEMBL4941', 'CHEMBL4208']
        ###################

        log(self.log_file, 'Getting chembl smiles files...')
        chembl_actives_smiles_files = get_files_with_suffix(self.chembl_smiles_path)  # load chembl smiles files
        log(self.log_file, 'Filtering chembl smiles...')

        self.filtered_chembl_active_files = []  # filter chembl smiles files only to those loaded from csv, also correlate path with chembl target
        for smile_file in chembl_actives_smiles_files:
            for target in self.targets_ids:
                if target in smile_file:
                    self.filtered_chembl_active_files.append([target, smile_file])
                    break
        self.chembl_data = {target: dict() for target in self.targets_ids}

        count = 0
        checkpoint = time()
        if not os.path.isfile('found_decoys_test/CHEMBL_values.pkl') or ignore_pickled is True:
            for target, smile_file in self.filtered_chembl_active_files:
                if count % 20 == 0:
                    log(self.log_file, f'Currently '
                        f'processed {round(count / len(self.filtered_chembl_active_files) * 100, 5)}% of targets')
                self.chembl_data[target] = self.load_chembl_targets(
                    smile_file)  # dict of targets as keys and list of compounds with smiles as values
                self.chembl_data[target] = self.add_descriptors(
                    self.chembl_data[target])  # add descriptors to compounds
                count += 1
            self.save_chembl_descriptors()
            log(self.log_file, f'Preparing targets took: {time() - checkpoint} s. Dumping them to csv.')
            self.dump_to_csv()
        else:
            log(self.log_file, f'Found chembl descriptors files, loading them...')
            self.chembl_data = pickle.load(open('found_decoys_test/CHEMBL_values.pkl', 'rb'))

        log(self.log_file, 'Getting ZINC files paths')
        self.zinc_database_smiles_files = sort_file_list_by_size(
            get_files_with_suffix(self.decoy_base_folder_path, suffix='.smi'))
        log(self.log_file, f'Starting parallelized job on 50 cores.')
        x = len([(file, dict({target: self.chembl_data[target]})) for target in self.chembl_data for file in
                 self.zinc_database_smiles_files])
        log(self.log_file, f'Number of tasks to do: {x}')
        Parallel(n_jobs=50)(delayed(ZINC_comparer)(file=file, chembl_data=dict({target: self.chembl_data[target]}),
                                                   threshold_dict=self.threshold_dicts,
                                                   max_decoys=self.max_decoys_for_active,
                                                   log_file=self.log_file) for target in self.chembl_data
                            for file in self.zinc_database_smiles_files)
        log(self.log_file, f'Whole pipeline took {time() - self.global_start_time} s.')

    def save_chembl_descriptors(self):
        pickle.dump(self.chembl_data, open('found_decoys_test/CHEMBL_values.pkl', 'wb'))

    def dump_to_csv(self):
        cols = ['chembl_target_id', 'ligand_id', 'HBD', 'HBA', 'rotates', 'weight', 'logp']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        i = 1
        for c_target in self.chembl_data:
            log(self.log_file, f'Currently on {i}/{len(self.chembl_data)}')
            i += 1
            for c_ligand in self.chembl_data[c_target]:
                r = [self.chembl_data[c_target][c_ligand][i] for i in ['HBD', 'HBA', 'rotates', 'weight', 'logp']]
                r = [c_target, c_ligand] + r
                output_csv.loc[index] = r
                index += 1
        output_csv.to_csv('found_decoys_test/CHEMBL_values.csv')

    def get_chembl_ids_from_csv(self, csv_path):
        master_table = pd.read_csv(csv_path, index_col=0)
        master_table = master_table[master_table['In_PDBBIND'] == 0]
        targets = master_table.index.tolist()
        return targets

    def load_chembl_targets(self, smiles_path):
        """
        Load whole smiles file as dict
        :param smiles_path: Path to smiles folder
        :return: Dict similar to those returned by load_smiles_in_batches
        """
        targets = dict()
        for smiles in self.load_smiles_in_batches(smiles_path, skip_lines=0):
            targets.update(smiles)  # list of tuples (smile, id)
        return targets


class ZINC_comparer(Decoy_finder):
    def __init__(self, file, chembl_data, threshold_dict, max_decoys, log_file):
        super().__init__()
        self.max_decoys_for_ligand = max_decoys
        self.log_file = log_file

        self.file_input = file
        self.file_name = file.split('/')[-1]
        log(self.log_file, f'Current ZINC file processed: {self.file_name}, target {list(chembl_data.keys())[0]}')
        # threshold_dict is a list of threshold_dicts
        self.results_zinc = [{} for _ in range(len(threshold_dict))]
        self.all_results = [{target: None for target in chembl_data} for _ in range(len(threshold_dict))]
        self.all_potential_decoys = [{target: None for target in chembl_data} for _ in range(len(threshold_dict))]
        self.found = [{target: {j: 0 for j in chembl_data[target]} for target in chembl_data} for _ in range(len(threshold_dict))]
        self.ignore_set = [set() for _ in range(len(threshold_dict))]

        self.count_observer = dict()
        for tar in chembl_data:
            self.count_observer[tar] = dict()
            for lig in tar:
                self.count_observer[tar][lig] = 0

        self.line_count = 0
        self.find_decoys(chembl_data, threshold_dict)

    def find_decoys(self, chembl_data, threshold_dict):
        checkpoint = time()
        with open(self.file_input, 'r') as handle:
            handle.readline()
            for line in handle:
                batch = dict()
                smile, smile_id = line.strip().split()
                batch[smile_id] = {'smile': smile, 'fingerprint': self.make_fingerprint(smile)}
                self.line_count += 1
                decoys_batch = batch
                decoys_batch = self.add_descriptors(decoys_batch)
                for threshold_idx in range(len(threshold_dict)):  # iterate over thresholds
                    for target in chembl_data:  # iterate over chembl targets
                        if target not in self.ignore_set[threshold_idx]:
                            data_from = {chembl_data[target][i] for i in set(chembl_data[target].keys())-self.ignore_set[idx]}  # remove ligands for which threshold was reached
                            # compare results in respect to desired thresholds and update
                            results, potential_decoys, found, chosen_ligand = self.compare_values(data_from, decoys_batch,
                                                                                   threshold_dict[threshold_idx],
                                                                                   self.all_results[threshold_idx][
                                                                                       target],
                                                                                   self.all_potential_decoys[
                                                                                       threshold_idx][target])
                            # update if there are any results
                            try:
                                if found > 0:
                                    try:
                                        self.results_zinc[threshold_idx].update(decoys_batch)
                                    except:
                                        self.results_zinc[threshold_idx] = decoys_batch
                                self.all_results[threshold_idx][target] = results
                                self.all_potential_decoys[threshold_idx][target] = potential_decoys

                                for idx in range(threshold_idx, len(threshold_dict)):   #  for each threshold
                                    self.found[threshold_idx][target][chosen_ligand] += found   #  add 1 to a ligand to which given ZINC matched
                                if self.found[threshold_idx][target][chosen_ligand] >= self.max_decoys_for_ligand:  # look if max decoys for ligands was found
                                    for idx in range(threshold_idx, len(threshold_dict)):   # in each threshold
                                        self.ignore_set[idx].add(chosen_ligand)
                                    log(self.log_file,
                                        f'Reached threshold {threshold_idx} for {target} in {self.file_name}')
                            except:
                                pass
                if self.line_count % 2500 == 0:
                    log(self.log_file, f'Currently processed {self.line_count} smiles for {self.file_name}.')
                b = True
                for idx in range(len(self.ignore_set)):
                    if len(self.ignore_set[idx]) < self.max_decoys_for_ligand:
                        b = False
                if b:
                    log(self.log_file, f'Found max decoy number for all targets in {self.file_name}')
                    break
        for threshold_idx in range(len(threshold_dict)):
            if sum(self.found[threshold_idx].values()) > 0:
                log(self.log_file, f'{self.found[threshold_idx]} found for {threshold_idx} in {self.file_name}')
                self.dump_to_csv(self.file_name, self.results_zinc[threshold_idx], threshold_idx)
                log(self.log_file, f'Saved potential decoys values for {self.file_name} threshold {threshold_idx}')
                self.dump_difference_to_csv(self.file_name, self.all_results[threshold_idx], threshold_idx)
                log(self.log_file,
                    f'Saved potential decoys and CHEMBL ligands differences for {self.file_name} threshold {threshold_idx}')
            else:
                s = f"No potential decoys found for {self.file_name} for threshold {str(threshold_idx)} \n"
                log(self.log_file, s)
        log(self.log_file,
            f'Finished search for {self.file_name} {list(chembl_data.keys())[0]} in {time() - checkpoint} s.')

    def dump_to_csv(self, filename, zinc_data, ix):
        cols = ['zinc', 'HBD', 'HBA', 'rotates', 'weight', 'logp']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        for c_zinc in zinc_data:
            r = [zinc_data[c_zinc][i] for i in ['HBD', 'HBA', 'rotates', 'weight', 'logp']]
            r = [c_zinc] + r
            output_csv.loc[index] = r
            index += 1
        if os.path.isfile(os.path.join('found_decoys_test', f'{filename}_thresholds_{ix}_values.csv')):
            other_csv = pd.read_csv(f'found_decoys_test/{filename}_thresholds_{ix}_values.csv', index_col=0)
            concated = pd.concat([other_csv, output_csv]).drop_duplicates().reset_index(drop=True)
            concated.to_csv(f'found_decoys_test/{filename}_thresholds_{ix}_values.csv')
        else:
            output_csv.to_csv(f'found_decoys_test/{filename}_thresholds_{ix}_values.csv')

    def dump_difference_to_csv(self, filename, difference_dict, ix):
        cols = ['target_id', 'ligand_id', 'zinc', 'HBD', 'HBA', 'rotates', 'weight', 'logp', 'Tanimoto',
                'murco_scaffold']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        for c_target in difference_dict:
            if type(difference_dict[c_target]) is dict:
                for c_ligand in difference_dict[c_target]:
                    for c_zinc in difference_dict[c_target][c_ligand]:
                        r = [difference_dict[c_target][c_ligand][c_zinc][i] for i in cols[3:]]
                        r = [c_target, c_ligand, c_zinc] + r
                        output_csv.loc[index] = r
                        index += 1

        if os.path.isfile(os.path.join('found_decoys_test', f'{filename}_thresholds_{ix}_differences.csv')):
            other_csv = pd.read_csv(f'found_decoys_test/{filename}_thresholds_{ix}_differences.csv', index_col=0)
            concated = pd.concat([other_csv, output_csv]).drop_duplicates().reset_index(drop=True)
            concated.to_csv(f'found_decoys_test/{filename}_thresholds_{ix}_differences.csv')
        else:
            output_csv.to_csv(f'found_decoys_test/{filename}_thresholds_{ix}_differences.csv')

    def compare_values(self, targets_from, targets_to, threshold_dict, result=None, potential_decoys=None):
        #  chcemy, aby monitorował po drodze, czy nie dobił do thresholdu znalezionych decoyidla danego ligandu
        #  no i jeśli decoy pasuje dla jednego ligandu to nie szukamy dalej
        if result is None:
            result = {t2: dict() for t2 in targets_from}
            potential_decoys = {t2: [] for t2 in targets_from}
        if threshold_dict is None:
            raise Exception('No cutoff values specified.')
        found = 0
        found_set = set()
        leave = 0
        for t1 in targets_to:
            target_to = targets_to[t1]
            for t2 in targets_from:
                if t1 != t2 and not leave:
                    target_from = targets_from[t2]
                    difference = dict()
                    difference['Tanimoto'] = DataStructs.FingerprintSimilarity(target_from['fingerprint'],
                                                                               target_to['fingerprint'])
                    difference['HBD'] = target_from['HBD'] - target_to['HBD']
                    difference['HBA'] = target_from['HBA'] - target_to['HBA']
                    difference['rotates'] = target_from['rotates'] - target_to['rotates']
                    difference['weight'] = target_from['weight'] - target_to['weight']
                    difference['logp'] = target_from['logp'] - target_to['logp']
                    difference['murco_scaffold'] = DataStructs.FingerprintSimilarity(target_from['murco_scaffold'],
                                                                                     target_to['murco_scaffold'])
                    merge = True
                    for i in difference:
                        dif = threshold_dict[i] - abs(difference[i])
                        if i == 'Tanimoto':
                            dif = difference[i] - threshold_dict[i]
                            if dif < 0:
                                merge = False
                                break
                        elif i in ['weight', 'logp']:
                            if difference[i] > target_from[i] * threshold_dict[i] or difference[i] < -(
                                    target_from[i] * threshold_dict[i]):
                                merge = False
                                break
                        else:
                            if dif < 0:
                                merge = False
                                break
                    if merge:
                        result[t2][t1] = difference
                        potential_decoys[t2].append((target_to['smile'], t1))
                        if t2 not in found_set:
                            found += 1
                            found_set.add(t2)
                            leave = 1
                            break
                        chosen_ligand = t2
            if leave:
                break
        return result, potential_decoys, found, chosen_ligand

    def save_decoys_to_file(self, file_name, all_potential_decoys, ix):
        with open(f'found_decoys_test/{file_name}-thresholds_{ix}-potential_decoys.txt', 'w') as handle:
            for target in all_potential_decoys:
                handle.write(f'>>>{target}\n')
                seen = set()
                for target_ligand in all_potential_decoys[target]:
                    for decoy in all_potential_decoys[target][target_ligand]:
                        if decoy[1] not in seen:
                            seen.add(decoy[1])
                            handle.write(f'{decoy[1]}\n')


def chunk(it, size):
    it = iter(it)
    return iter(lambda: tuple(islice(it, size)), ())


def sort_file_list_by_size(file_list):
    pairs = []
    for file in file_list:
        # Get size and add to list of tuples.
        size = os.path.getsize(file)
        pairs.append((size, file))

    # Sort list of tuples by the first element, size.
    pairs.sort(key=lambda s: s[0])

    # Display pairs.
    result = [pair[1] for pair in pairs]
    return result


parser = argparse.ArgumentParser()
parser.add_argument('--zinc', default='raw_data/ZINC_test', help="Folder with ZINC files")
parser.add_argument('--chembl_csv', default='actives_number_sampling/targets_after_fingerprint_similarity0_tc0.95.csv',
                    help="CHEMBL csv")
parser.add_argument('--chembl_smiles', default='chembl_smiles', help="Folder with CHEMBL ligands smiles")
parser.add_argument('--max_decoys', default=200, help="Max decoys gor target per ZINC file")
parser.add_argument('--log_file', default='logs_test.txt', help="Log to file?")
parser.add_argument('--ignore_pickled', type=lambda x: bool(distutils.util.strtobool(x)), default=True,
                    help="Ignore pickled CHEMBL values save if found in dest folder?")

if __name__ == '__main__':
    args = parser.parse_args()
    chembl = args.chembl_csv
    chembl_smiles = args.chembl_smiles
    decoy_path = args.zinc

    # max difference between target and new decoy,
    # except TC which indicates min. value of that descriptor that
    thresholds = [
        {'Tanimoto': 0.85, 'HBD': 0, 'HBA': 0, 'rotates': 0, 'weight': 0.05, 'logp': 0.05, 'murco_scaffold': 0.7},
        {'Tanimoto': 0.8, 'HBD': 2, 'HBA': 2, 'rotates': 2, 'weight': 0.05, 'logp': 0.05, 'murco_scaffold': 0.7},
        {'Tanimoto': 0.8, 'HBD': 3, 'HBA': 3, 'rotates': 3, 'weight': 0.08, 'logp': 0.08, 'murco_scaffold': 0.7},
        {'Tanimoto': 0.8, 'HBD': 4, 'HBA': 4, 'rotates': 4, 'weight': 0.10, 'logp': 0.10, 'murco_scaffold': 0.7},
        {'Tanimoto': 0.75, 'HBD': 5, 'HBA': 5, 'rotates': 4, 'weight': 0.15, 'logp': 0.15, 'murco_scaffold': 0.7}
    ]
    decoys_per_target_per_file = args.max_decoys

    Chembl_loader(chembl, chembl_smiles, decoy_path, thresholds, decoys_per_target_per_file, args.log_file,
                  args.ignore_pickled)
