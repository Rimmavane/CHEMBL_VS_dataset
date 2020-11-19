from os.path import isfile, basename
import pandas as pd
import argparse
import pickle
import random
from utils import log, get_files_with_suffix, make_fingerprint
from paths_and_settings import *

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, MACCSkeys
import distutils.util
from rdkit import DataStructs
import rdkit.Chem.Scaffolds.MurckoScaffold as murco


class Decoy_finder:
    def __init__(self):
        pass

    def add_descriptors(self, targets_dict):
        for target in targets_dict:
            desc = self.make_descriptors_for_smile(targets_dict[target]['smile'])
            targets_dict[target].update(desc)
        return targets_dict

    @staticmethod
    def load_smiles(file_path, skip_lines=0):
        with open(file_path, 'r') as handle:
            for _ in range(skip_lines):
                next(handle)
            smiles = dict()
            i = 0
            for line in handle:
                smile, smile_id = line.strip().split()
                smiles[smile_id] = {'smile': smile, 'fingerprint': make_fingerprint(smile)}
                i += 1
        return smiles

    @staticmethod
    def make_descriptors_for_smile(smile, fp='fp2'):
        """
        Make physical descriptors:
        - molecular weight,
        - the number of rotational bonds,
        - total hydrogen bond donors (HBDs),
        - total hydrogen bond acceptors (HBAs),
        - octanol–water partition coefficient (log P)
        :param smile: SMILE representation of compound
        :aram fp: Fingerprint type to use, 'fp2', 'morg2' and 'maccs' are availible
        :return: Dictionary with physical descriptors of smile
        """
        descriptors = dict()
        x = Chem.MolFromSmiles(smile)
        descriptors['HBD'] = rdMolDescriptors.CalcNumHBD(x)
        descriptors['HBA'] = rdMolDescriptors.CalcNumLipinskiHBA(x)
        descriptors['rotates'] = rdMolDescriptors.CalcNumRotatableBonds(x)
        descriptors['weight'] = rdMolDescriptors.CalcExactMolWt(x)
        descriptors['logp'] = rdMolDescriptors.CalcCrippenDescriptors(x)[0]
        if fp == 'fp2':
            descriptors['murco_scaffold'] = Chem.RDKFingerprint(murco.GetScaffoldForMol(x))
        elif fp == 'morg2':
            descriptors['murco_scaffold'] = AllChem.GetMorganFingerprintAsBitVect(murco.GetScaffoldForMol(x), 2, useFeatures=True)
        else:
            descriptors['murco_scaffold'] = MACCSkeys.GenMACCSKeys(murco.GetScaffoldForMol(x))
        return descriptors

    @staticmethod
    def parse_dude_ism(file_path):
        with open(file_path, 'r') as handle:
            smiles = dict()
            for line in handle:
                l = line.strip().split()
                smiles[l[-1]] = {'smile': l[0], 'fingerprint': make_fingerprint(l[0])}
        return smiles


class Chembl_loader(Decoy_finder):
    def __init__(self, chembl_id, chembl_smi_path, decoy_folder_path, th_dict, max_decoys, output_folder='', ignore_pickled=True, randomize_ligands=True):
        if output_folder == '':
            self.output_folder = join(PROJECT_HOME, 'found_decoys')
        else:
            self.output_folder = output_folder
        super().__init__()
        self.global_start_time = time.time()
        self.chembl_smiles_path = str(chembl_smi_path)
        self.decoy_base_folder_path = str(decoy_folder_path)
        self.threshold_dict = th_dict
        self.max_decoys_for_ligand = int(max_decoys)
        self.chembl_id = str(chembl_id)
        self.target_ids = [self.chembl_id]

        self.filtered_chembl_active_files = [[target, join(self.chembl_smiles_path, f'{target}_filtered_active.smi')]
                                             for target in self.target_ids]
        self.chembl_data = {target: dict() for target in self.target_ids}

        create_folder_if_not_existent(self.output_folder)

        checkpoint = time.time()
        if not isfile(join(self.output_folder, f'{self.chembl_id}_descriptors_values.pkl')) or ignore_pickled is True:
            for target, smile_file in self.filtered_chembl_active_files:
                self.chembl_data[target] = self.load_chembl_ligands(smile_file)  # dict of targets as keys and list of compounds with smiles as values
                self.chembl_data[target] = self.add_descriptors(self.chembl_data[target])  # add descriptors to compounds
            self.save_chembl_descriptors()
            log(f'Preparing targets took: {time.time() - checkpoint} s. Dumping them to csv.')
            self.dump_to_csv()
        else:
            log(f'Found chembl descriptors files, loading them...')
            self.chembl_data = pickle.load(open(join(self.output_folder, f'{self.chembl_id}_descriptors_values.pkl'), 'rb'))

        log('Getting ZINC files paths')
        self.zinc_database_smiles_files = get_files_with_suffix(self.decoy_base_folder_path, suffix='.smi')
        x = len([(file, dict({target: self.chembl_data[target]})) for target in self.chembl_data
                 for file in self.zinc_database_smiles_files])
        log(f'Number of tasks to do: {x}')
        for file in self.zinc_database_smiles_files:
            for target in self.chembl_data:
                if randomize_ligands:  # randomize ligand order so if decoys are searched in many files the decoys found, which suffice thresholds for more than one active, will be spread more evenly
                    random_ligands = list(self.chembl_data[target].items())
                    random.shuffle(random_ligands)
                    self.chembl_data[target] = dict(random_ligands)
                ZINC_comparer(chembl_id=target,
                              decoys_file_path=file,
                              chembl_data=dict({target: self.chembl_data[target]}),
                              threshold_dict=self.threshold_dict,
                              max_decoys=self.max_decoys_for_ligand,
                              output_folder=self.output_folder)
        log(f'Whole pipeline for {self.chembl_id} took {time.time() - self.global_start_time} s.')

    def save_chembl_descriptors(self):
        pickle.dump(self.chembl_data, open(join(self.output_folder, f'{self.chembl_id}_descriptors_values.pkl'), 'wb'))

    def dump_to_csv(self):
        cols = ['chembl_target_id', 'ligand_id', 'HBD', 'HBA', 'rotates', 'weight', 'logp']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        i = 1
        for c_target in self.chembl_data:
            log(f'Currently on {i}/{len(self.chembl_data)}')
            i += 1
            for c_ligand in self.chembl_data[c_target]:
                r = [self.chembl_data[c_target][c_ligand][i] for i in ['HBD', 'HBA', 'rotates', 'weight', 'logp']]
                r = [c_target, c_ligand] + r
                output_csv.loc[index] = r
                index += 1
        output_csv.to_csv(join(self.output_folder, f'{self.chembl_id}_descriptors_values.csv'))

    def get_chembl_ids_from_csv(self, csv_path):
        master_table = pd.read_csv(csv_path, index_col=0)
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
    def __init__(self, chembl_id, decoys_file_path, chembl_data, threshold_dict, max_decoys, output_folder):
        super().__init__()
        self.checkpoint = time.time()
        self.max_decoys_for_ligand = max_decoys
        self.output_folder = output_folder

        self.chembl_id = chembl_id
        self.zinc_file_path = decoys_file_path   # path to decoy files
        self.file_name = basename(decoys_file_path)  # extract filename from path
        log(f'Current ZINC file processed: {self.file_name}, target {list(chembl_data.keys())[0]}')
        # self.results_zinc = [dict() for _ in range(len(threshold_dict))]
        self.all_results = {target: {ligand: dict() for ligand in chembl_data[target]} for target in chembl_data}
        self.found = {target: {ligand: 0 for ligand in chembl_data[target]} for target in chembl_data}
        self.ignore_targets = set()
        self.ignore_set_chembl = set()
        self.ignore_set_zinc = {target: set() for target in chembl_data}
        self.threshold_counts = 0

        self.line_count = 0
        assert len(chembl_data) == 1
        self.find_decoys(chembl_data, threshold_dict)

    def find_decoys(self, chembl_data, threshold_dict):
        with open(self.zinc_file_path, 'r') as handle:
            handle.readline()
            for line in handle:  # read smile by smile
                self.line_count += 1
                decoys_batch = self.read_smile_fp(line)
                if len(self.ignore_targets) == len(chembl_data):
                    log(f'Reached threshold for all ligands in {self.chembl_id} target')
                    break
                for target in chembl_data:  # chembl_data[target] == {LIGAND_ID1:{descriptors}, LIGAND_ID2: {....}...}
                    if target not in self.ignore_targets:
                        for ligand in chembl_data[target]:  # get into CHMEBL_ID dict --> iterate over ligands
                            if len(self.ignore_set_chembl) == len(chembl_data[target]):
                                self.ignore_targets.add(target)
                                break
                            if ligand not in self.ignore_set_chembl and self.found[target][ligand] < self.max_decoys_for_ligand:
                                difference, good = self.compare_values(ligand, chembl_data[target], decoys_batch, threshold_dict)  # ligand_id, chembl_id_dict, zinc_smiles_dict, thresholds
                                # upper function returns:
                                # dict of differences between each ligand and zinc
                                # update if there are any results
                                if good == 1 and difference['zinc_id'] not in self.ignore_set_zinc[target]:
                                    self.all_results[target][ligand][difference['zinc_id']] = difference
                                    self.ignore_set_zinc[target].add(difference['zinc_id'])
                                    self.found[target][ligand] += 1
                                    if self.found[target][ligand] == self.max_decoys_for_ligand:
                                        self.ignore_set_chembl.add(ligand)
                                        log(f'Reached found decoys threshold for ligand {ligand} in {target} in {self.file_name} {len(self.ignore_set_chembl)}/{len(chembl_data[target])}')
                                    self.threshold_counts += 1
                                    break
                if self.line_count % 2500 == 0:
                    log(f'Currently processed {self.line_count} smiles from {self.file_name}.')
        if self.threshold_counts > 0:
            log(f'{self.threshold_counts} found for {self.chembl_id} in {self.file_name}')
            self.dump_difference_to_csv(self.file_name, self.all_results)
            log(f'Saved potential decoys and CHEMBL ligands differences for {self.file_name}')
        else:
            message = f"No potential decoys found for {self.file_name}"
            log(message)
        log(f'Finished search for {self.file_name} {list(chembl_data.keys())[0]} in {time.time() - self.checkpoint} s.')

    def read_smile_fp(self, line):
        decoys_batch = dict()
        smile, smile_id = line.strip().split()
        decoys_batch[smile_id] = {'smile': smile, 'fingerprint': make_fingerprint(smile)}
        decoys_batch = self.add_descriptors(decoys_batch)
        return decoys_batch

    def dump_difference_to_csv(self, filename, difference_dict):
        cols = ['target_id', 'ligand_id', 'zinc', 'ligand_smile', 'zinc_smile', 'HBD', 'HBA',
                'rotates', 'weight', 'logp', 'murco_scaffold']
        output_csv = pd.DataFrame(columns=cols)
        index = 0
        for c_target in difference_dict:
            for c_ligand in difference_dict[c_target]:
                for c_zinc in difference_dict[c_target][c_ligand]:
                    r = [difference_dict[c_target][c_ligand][c_zinc][i] for i in cols[3:]]
                    r = [c_target, c_ligand, c_zinc] + r
                    output_csv.loc[index] = r
                    index += 1
        if isfile(join(self.output_folder, f'{filename}_differences.csv')):
            other_csv = pd.read_csv(join(self.output_folder, f'{filename}_differences.csv'), index_col=0)
            concated = pd.concat([other_csv, output_csv]).drop_duplicates().reset_index(drop=True)
            concated.to_csv(join(self.output_folder, f'{filename}_differences.csv'))
        else:
            output_csv.to_csv(join(self.output_folder, f'{filename}_differences.csv'))

    @staticmethod
    def compare_values(ligand_id, chembl_data_target, decoys_batch, threshold_dict):
        if threshold_dict is None:
            raise Exception('No cutoff values specified.')
        difference = dict()
        good = 0
        for zinc_id in decoys_batch:  # get that one ZINC ID
            difference = dict()  # calculate differences
            good = 1
            target_to = decoys_batch[zinc_id]  # ZINC ID descriptors values
            target_from = chembl_data_target[ligand_id]   # LIGAND ID descritors values
            difference['HBD']            = target_from['HBD'] - target_to['HBD']
            difference['HBA']            = target_from['HBA'] - target_to['HBA']
            difference['rotates']        = target_from['rotates'] - target_to['rotates']
            difference['weight']         = (target_from['weight'] - target_to['weight'])/target_from['weight']
            difference['logp']           = (target_from['logp'] - target_to['logp'])/(0.001 + target_from['logp'])
            difference['murco_scaffold'] = DataStructs.FingerprintSimilarity(target_from['murco_scaffold'],
                                                                             target_to['murco_scaffold'])
            for i in difference:    # compare differences with thresholds
                if i in ['weight']:
                    if abs(difference[i]) > threshold_dict[i]:
                        good = 0
                        break
                elif i in ['logp']:
                    if SOFT_LOGP_THRESHOLDS:
                        if abs(difference[i])-(1-(abs(target_from['logp'])/SOFT_LOGP_PARAMETER)) > threshold_dict[i]:   # adaptive comparison
                            good = 0
                            break
                    else:
                        if abs(difference[i]) > threshold_dict[i]:  # hard comaprison
                            good = 0
                            break
                else:
                    dif = threshold_dict[i] - abs(difference[i])  # HBD, HBA, rotates, murco_scaffold
                    if dif < 0:
                        good = 0
                        break
            if good:
                difference['ligand_id'] = ligand_id
                difference['ligand_smile'] = target_from['smile']
                difference['zinc_id'] = zinc_id
                difference['zinc_smile'] = target_to['smile']
            break
        return difference, good


parser = argparse.ArgumentParser()
parser.add_argument('--decoys_path', default='raw_data/ZINC', help="Folder with ZINC files")
parser.add_argument('--chembl_id', default='', help="Target for which the decoy search will be executed.")
parser.add_argument('--chembl_smiles', default=CHEMBL_SMILES_FOLDER, help="Folder with CHEMBL ligands smiles. By default same as CHEMBL_SMILES_FOLDER parameter in path_and_settings script")
parser.add_argument('--max_decoys', default=100, help="Max decoys for each ligand for a target per ZINC file")
parser.add_argument('--output_folder', default='', help="Where to save found ")
parser.add_argument('--ignore_pickled', type=lambda x: bool(distutils.util.strtobool(x)), default=True, help="Ignore pickled CHEMBL values  if found in destination folder and save anyway")
parser.add_argument('--randomize_ligands', type=lambda x: bool(distutils.util.strtobool(x)), default=True, help="Randomize ligands order in actives set")

if __name__ == '__main__':
    args = parser.parse_args()
    assert args.chembl_id

    Chembl_loader(chembl_id=args.chembl_id, chembl_smi_path=args.chembl_smiles, decoy_folder_path=args.decoys_path, th_dict=THRESHOLDS_DICT,
                  max_decoys=args.max_decoys, output_folder=args.output_folder, ignore_pickled=args.ignore_pickled, randomize_ligands=args.randomize_ligands)
