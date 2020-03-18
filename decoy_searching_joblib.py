import os
import pandas as pd
from time import time
import argparse
import pickle
import distutils.util

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from rdkit import DataStructs
import rdkit.Chem.Scaffolds.MurckoScaffold as murco
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

    def load_smiles(self, file_path, skip_lines=0):
        with open(file_path, 'r') as handle:
            for _ in range(skip_lines):
                next(handle)
            smiles = dict()
            i = 0
            for line in handle:
                smile, smile_id = line.strip().split()
                smiles[smile_id] = {'smile': smile, 'fingerprint': self.make_fingerprint(smile)}
                i += 1
            return smiles

    def make_fingerprint(self, smile):
        """
        Make ECFP4 fingerprint for given smile
        :param smile: Compound smile
        :return: ECFP4 fingerprint
        """
        try:
            x = Chem.MolFromSmiles(smile)
            x = AllChem.GetMorganFingerprint(x, 2)
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
        descriptors['murco_scaffold'] = AllChem.GetMorganFingerprint(murco.GetScaffoldForMol(x), 2)
        return descriptors


class Chembl_loader(Decoy_finder):
    def __init__(self, chembl_csv, chembl_smi_path, decoy_folder_path, th_dicts,
                 max_decoys, log_file='', ignore_pickled=True):
        super().__init__()
        self.log_file = log_file
        self.n_jobs = os.cpu_count()-1
        self.global_start_time = time()
        self.chembl_smiles_path = chembl_smi_path
        self.decoy_base_folder_path = decoy_folder_path
        self.threshold_dicts = th_dicts
        self.max_decoys_for_ligand = max_decoys
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'): pass
        log(self.log_file, f"Availible cores: {self.n_jobs}")

        log(self.log_file, 'Loading chembl ids...')
        if type(chembl_csv) != list:
            self.targets_ids = self.get_chembl_ids_from_csv(chembl_csv)  # load chembl ids
        else:
            self.targets_ids = chembl_csv

        ### TEST SUBJECT ##
        # self.targets_ids = ['CHEMBL4941', 'CHEMBL4208']
        self.targets_ids = ['CHEMBL3784']
        ###################

        log(self.log_file, 'Getting chembl smiles files...')
        chembl_actives_smiles_files = get_files_with_suffix(self.chembl_smiles_path, suffix='_active.smi')  # load chembl smiles files
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
        if not os.path.isfile('found_decoys/CHEMBL_values.pkl') or ignore_pickled is True:
            for target, smile_file in self.filtered_chembl_active_files:
                if count % 20 == 0:
                    log(self.log_file, f'Currently processed {round(count / len(self.filtered_chembl_active_files) * 100, 5)}% of targets')
                self.chembl_data[target] = self.load_chembl_ligands(
                    smile_file)  # dict of targets as keys and list of compounds with smiles as values
                self.chembl_data[target] = self.add_descriptors(
                    self.chembl_data[target])  # add descriptors to compounds
                count += 1
            self.save_chembl_descriptors()
            log(self.log_file, f'Preparing targets took: {time() - checkpoint} s. Dumping them to csv.')
            self.dump_to_csv()
        else:
            log(self.log_file, f'Found chembl descriptors files, loading them...')
            self.chembl_data = pickle.load(open('found_decoys/CHEMBL_values.pkl', 'rb'))

        log(self.log_file, 'Getting ZINC files paths')
        self.zinc_database_smiles_files = get_files_with_suffix(self.decoy_base_folder_path, suffix='.smi')
        log(self.log_file, f'Starteing parallelized job on {2*os.cpu_count()} cores')
        x = len([(file, dict({target: self.chembl_data[target]})) for target in self.chembl_data
                 for file in self.zinc_database_smiles_files])
        log(self.log_file, f'Number of tasks to do: {x}')
        # Parallel(n_jobs=2*os.cpu_count())(delayed(ZINC_comparer)(file=file, chembl_data=dict({target: self.chembl_data[target]}),  # execute for each chembl target and zinc file
        #                                            threshold_dict=self.threshold_dicts,
        #                                            max_decoys=self.max_decoys_for_ligand,
        #                                            log_file=self.log_file) for target in self.chembl_data
        #                                                                   for file in self.zinc_database_smiles_files)
        Parallel(n_jobs=4)(delayed(ZINC_comparer)(file=file, chembl_data=dict({target: self.chembl_data[target]}),  # execute for each chembl target and zinc file
                                           threshold_dict=self.threshold_dicts,
                                           max_decoys=self.max_decoys_for_ligand,
                                           log_file=self.log_file) for target in self.chembl_data
                                                                  for file in self.zinc_database_smiles_files)

        log(self.log_file, f'Whole pipeline took {time() - self.global_start_time} s.')

    def save_chembl_descriptors(self):
        pickle.dump(self.chembl_data, open('found_decoys/CHEMBL_values.pkl', 'wb'))

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
        output_csv.to_csv('found_decoys/CHEMBL_values.csv')

    def get_chembl_ids_from_csv(self, csv_path):
        master_table = pd.read_csv(csv_path, index_col=0)
        master_table = master_table[master_table['In_PDBBIND'] == 0]
        targets = master_table.index.tolist()
        return targets

    def load_chembl_ligands(self, smiles_path):
        """
        Load whole smiles file as dict
        :param smiles_path: Path to smiles folder
        :return: Dict similar to those returned by load_smiles_in_batches
        """
        targets = self.load_smiles(smiles_path, skip_lines=0)  # list of tuples (smile, id)
        return targets


class ZINC_comparer(Decoy_finder):
    def __init__(self, file, chembl_data, threshold_dict, max_decoys, log_file):
        super().__init__()
        self.max_decoys_for_ligand = max_decoys
        self.log_file = log_file

        self.zinc_file_path = file   # path to zinc file
        self.file_name = os.path.basename(file) # extract filename from path
        log(self.log_file, f'Current ZINC file processed: {self.file_name}, target {list(chembl_data.keys())[0]}')
        #threshold_dict is a list of threshold_dicts
        # self.results_zinc = [dict() for _ in range(len(threshold_dict))]
        self.all_results = [{target: {ligand: dict() for ligand in chembl_data[target]} for target in chembl_data} for _ in range(len(threshold_dict))]
        self.found = [{target: {ligand: 0 for ligand in chembl_data[target]} for target in chembl_data} for _ in range(len(threshold_dict))]
        self.ignore_set_chembl = [set() for _ in range(len(threshold_dict))]
        self.ignore_set_zinc = {target: set() for target in chembl_data}
        self.threshold_counts = [0 for _ in range(len(threshold_dict))]

        self.line_count = 0
        self.find_decoys(chembl_data, threshold_dict)

    def find_decoys(self, chembl_data, threshold_dict):
        self.checkpoint = time()
        with open(self.zinc_file_path, 'r') as handle:
            handle.readline()
            line_counter = 0
            for line in handle:  # read smile by smile
                bail = False
                self.line_count += 1
                decoys_batch = dict()
                smile, smile_id = line.strip().split()
                decoys_batch[smile_id] = {'smile': smile, 'fingerprint': self.make_fingerprint(smile)}
                decoys_batch = self.add_descriptors(decoys_batch)
                for target in chembl_data:  # chembl_data[target] == {LIGAND_ID1:{descriptors}, Ligand_id2: {....}...}
                    for threshold_idx in range(len(threshold_dict)): #iterate over thresholds
                        for ligand in chembl_data[target]:  # get into CHMEBL_ID dict --> iterate over ligands
                            if ligand not in self.ignore_set_chembl[threshold_idx]:
                                #compare results in respect to desired thresholds and update
                                difference, good = self.compare_values(ligand, chembl_data[target], decoys_batch,   #ligand_id, chembl_id_dict, zinc_smiles_dict
                                                            threshold_dict[threshold_idx])   #threshold values)
                                # upper function returns:
                                # dict of differences between each ligand and zinc
                                # update if there are any results
                                if good == 1 and difference['zinc_id'] not in self.ignore_set_zinc[target]:    #if there was a match not seen before
                                    self.all_results[threshold_idx][target][ligand][difference['zinc_id']] = difference
                                    self.ignore_set_zinc[target].add(difference['zinc_id'])
                                    # self.results_zinc[threshold_idx].update(decoys_batch)   # mark ZINC as valid decoy
                                    # self.all_results[threshold_idx][target][ligand] = results #happens inside compare_values
                                    for idx in range(threshold_idx, len(threshold_dict)):
                                        self.found[idx][target][ligand] += 1
                                    if self.found[threshold_idx][target][ligand] >= self.max_decoys_for_ligand:
                                        for idx in range(threshold_idx, len(threshold_dict)):
                                            self.ignore_set_chembl[idx].add(ligand)
                                        log(self.log_file, f'Reached threshold {threshold_idx} for {target} for {ligand} in {self.file_name}')
                                    self.threshold_counts[threshold_idx] += 1
                                    bail = True
                                    break
                            if bail:
                                break
                    line_counter += 1
                    if line_counter % 2500 == 0:
                        log(self.log_file, f'Currently processed {line_counter} smiles for {self.file_name}.')
                    if bail:
                        break
        log(self.log_file, f'{self.threshold_counts}')
        for threshold_idx in range(len(threshold_dict)):
            if self.threshold_counts[threshold_idx] > 0:
                log(self.log_file, f'{self.threshold_counts[threshold_idx]} found for {threshold_idx} in {self.file_name}')
                self.dump_difference_to_csv(self.file_name, self.all_results[threshold_idx], threshold_idx)
                log(self.log_file, f'Saved potential decoys and CHEMBL ligands differences for {self.file_name} threshold {threshold_idx}')
            else:
                message = f"No potential decoys found for {self.file_name} for threshold {str(threshold_idx)}"
                log(self.log_file, message)
        log(self.log_file, f'Finished search for {self.file_name} {list(chembl_data.keys())[0]} in {time() - self.checkpoint} s.')

    def dump_difference_to_csv(self, filename, difference_dict, ix):
        cols = ['target_id', 'ligand_id', 'zinc', 'ligand_smile', 'zinc_smile', 'HBD', 'HBA',
                'rotates', 'weight', 'logp', 'Tanimoto', 'murco_scaffold']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        for c_target in difference_dict:
            for c_ligand in difference_dict[c_target]:
                for c_zinc in difference_dict[c_target][c_ligand]:
                    r = [difference_dict[c_target][c_ligand][c_zinc][i] for i in cols[3:]]
                    r = [c_target, c_ligand, c_zinc] + r
                    output_csv.loc[index] = r
                    index += 1
        if os.path.isfile(os.path.join('found_decoys', f'{filename}_thresholds_{ix}_differences.csv')):
            other_csv = pd.read_csv(f'found_decoys/{filename}_thresholds_{ix}_differences.csv', index_col=0)
            concated = pd.concat([other_csv, output_csv]).drop_duplicates().reset_index(drop=True)
            concated.to_csv(f'found_decoys/{filename}_thresholds_{ix}_differences.csv')
            # self.test_tanimoto_value(f'found_decoys/{filename}_thresholds_{ix}_differences.csv')
        else:
            output_csv.to_csv(f'found_decoys/{filename}_thresholds_{ix}_differences.csv')
            # self.test_tanimoto_value(f'found_decoys/{filename}_thresholds_{ix}_differences.csv')

    def compare_values(self, ligand_id, chembl_data_target, decoys_batch, threshold_dict):
        if threshold_dict is None:
            raise Exception('No cutoff values specified.')
        difference = dict()
        good = 0
        for zinc_id in decoys_batch:  # get that one ZINC ID
            difference = dict()  # calculate differences
            good = 1
            target_to = decoys_batch[zinc_id]  # ZINC ID descriptors values
            target_from = chembl_data_target[ligand_id]   # LIGAND ID descritors values
            difference['Tanimoto'] = DataStructs.TanimotoSimilarity(target_from['fingerprint'],
                                                                    target_to['fingerprint'])
            difference['HBD'] = target_from['HBD'] - target_to['HBD']
            difference['HBA'] = target_from['HBA'] - target_to['HBA']
            difference['rotates'] = target_from['rotates'] - target_to['rotates']
            difference['weight'] = (target_from['weight'] - target_to['weight'])/target_from['weight']
            difference['logp'] = (target_from['logp'] - target_to['logp'])/target_from['logp']
            difference['murco_scaffold'] = DataStructs.TanimotoSimilarity(target_from['murco_scaffold'],
                                                                          target_to['murco_scaffold'])
            for i in difference:    # compare differences with thresholds
                if i == 'Tanimoto':
                    dif = difference[i] - threshold_dict[i]
                    if dif <= 0 or difference[i] == 1:
                        good = 0
                        break
                elif i in ['weight', 'logp']:
                    if abs(difference[i]) >= threshold_dict[i]:
                        good = 0
                        break
                else:
                    dif = threshold_dict[i] - abs(difference[i])  # HBD, HBA, rotates, murco_scaffold
                    if dif <= 0:
                        good = 0
                        break
            if good:
                difference['ligand_id'] = ligand_id
                difference['ligand_smile'] = target_from['smile']
                difference['zinc_id'] = zinc_id
                difference['zinc_smile'] = target_to['smile']
            break
        return difference, good

    def test_tanimoto_value(self, f):
        c = pd.read_csv(f, index_col=0)
        x = c.loc[0, :]
        s2 = x['zinc_smile']
        s1 = x['ligand_smile']
        print(s1, s2)
        s1 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(s1), 2)
        s2 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(s2), 2)
        print(DataStructs.TanimotoSimilarity(s1, s2), x['Tanimoto'])
        assert round(DataStructs.TanimotoSimilarity(s1, s2), 5) == round(x['Tanimoto'], 5)

def sort_file_list_by_size(file_list):
    pairs = []
    for file in file_list:

        # Get size and add to list of tuples.
        size = os.path.getsize(file)
        pairs.append((size, file))

    # Sort list of tuples by the first element, size.
    pairs.sort(key=lambda s: s[0], reverse=True)

    # Display pairs.
    result = [pair[1] for pair in pairs]
    return result


parser = argparse.ArgumentParser()
parser.add_argument('--zinc', default='raw_data/ZINC_shuffled', help="Folder with ZINC files")
parser.add_argument('--chembl_csv', default='actives_number_sampling/targets_after_fingerprint_similarity5_tc0.95.csv', help="CHEMBL csv")
parser.add_argument('--chembl_smiles', default='chembl_smiles', help="Folder with CHEMBL ligands smiles")
parser.add_argument('--max_decoys', default=25, help="Max decoys gor target per ZINC file")
parser.add_argument('--log_file', default='logs_decoys.txt', help="Log to file?")
parser.add_argument('--ignore_pickled', type=lambda x: bool(distutils.util.strtobool(x)), default=True, help="Ignore pickled CHEMBL values save if found in dest folder?")

if __name__ == '__main__':
    args = parser.parse_args()
    chembl = args.chembl_csv
    chembl_smiles = args.chembl_smiles
    decoy_path = args.zinc

    # max difference between target and new decoy,
    # except TC which indicates min. value of that descriptor that
    # thresholds = [
    #     {'Tanimoto': 0.85, 'HBD': 0, 'HBA': 0, 'rotates': 0, 'weight': 0.05, 'logp': 0.05, 'murco_scaffold': 0.7},
    #     {'Tanimoto': 0.8, 'HBD': 2, 'HBA': 2, 'rotates': 2, 'weight': 0.05, 'logp': 0.05, 'murco_scaffold': 0.7},
    #     {'Tanimoto': 0.8, 'HBD': 4, 'HBA': 4, 'rotates': 2, 'weight': 0.08, 'logp': 0.08, 'murco_scaffold': 0.7},
    #     {'Tanimoto': 0.8, 'HBD': 4, 'HBA': 4, 'rotates': 4, 'weight': 0.15, 'logp': 0.12, 'murco_scaffold': 0.7},
    #     {'Tanimoto': 0.75, 'HBD': 5, 'HBA': 5, 'rotates': 4, 'weight': 0.20, 'logp': 0.20, 'murco_scaffold': 0.7}s,
    #     {'Tanimoto': 0.5, 'HBD': 6, 'HBA': 6, 'rotates': 6, 'weight': 0.25, 'logp': 0.25, 'murco_scaffold': 0.75}
    # ]
    thresholds = [
        {'Tanimoto': 0.75, 'HBD': 1, 'HBA': 1, 'rotates': 1, 'weight': 0.15, 'logp': 0.15, 'murco_scaffold': 0.6},
        {'Tanimoto': 0.7, 'HBD': 2, 'HBA': 2, 'rotates': 2, 'weight': 0.15, 'logp': 0.15, 'murco_scaffold': 0.6},
        {'Tanimoto': 0.65, 'HBD': 3, 'HBA': 3, 'rotates': 3, 'weight': 0.15, 'logp': 0.15, 'murco_scaffold': 0.6}]
    decoys_per_target_per_file = args.max_decoys

    Chembl_loader(chembl_csv=chembl, chembl_smi_path=chembl_smiles, decoy_folder_path=decoy_path, th_dicts=thresholds,
                  max_decoys=decoys_per_target_per_file, log_file=args.log_file, ignore_pickled=args.ignore_pickled)
