import os
import pandas as pd
from tanimoto_coef_allvsall import log
import active_vs_inactive_bias as enrichment
from joblib import Parallel, delayed
from rdkit import Chem, DataStructs
import pickle
from rdkit.Chem import AllChem


def create_folder_if_not_existent(path):
    if not os.path.exists(path):
        os.makedirs(path)


def load_smiles(filepath):
    smiles_dict = dict()
    # if path to file provided
    try:
        with open(filepath, 'r') as smiles_file:
            for line in smiles_file:
                line = line.strip().split()
                smiles_dict[line[1]] = line[0]
    # if only chembl name provided
    except FileNotFoundError:
        file = \
            [i for i in os.listdir('chembl_smiles') if
             os.path.isfile(os.path.join('chembl_smiles', i)) and filepath in i][0]
        with open(f'chembl_smiles/{file}', 'r') as smiles_file:
            for line in smiles_file:
                line = line.strip().split()
                smiles_dict[line[1]] = line[0]
    return smiles_dict


# prepare decoys csv with :
# tail -n+2 * -q > decoys_amalgam.csv
class Prepare:
    def __init__(self, path_found_decoys='found_decoys', top=50,
                 storage='decoys_matrices', log_file=''):
        self.log_file = log_file
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'):
                pass

        self.storage = storage
        create_folder_if_not_existent(storage)

        log(log_file, 'Beginning loading found decoys...')
        self.found_decoys = pd.read_csv(os.path.join(path_found_decoys, 'decoys_amalgam.csv'),
                                        header=None,
                                        names=['', 'target_id', 'ligand_id', 'zinc', 'zinc_smile', 'HBD', 'HBA',
                                               'rotates', 'weight', 'logp', 'Tanimoto', 'murco_scaffold'])
        self.found_decoys = self.found_decoys.drop_duplicates(subset="zinc_smile")
        self.found_decoys = self.found_decoys.sort_values(['Tanimoto', 'HBD', 'HBA', 'rotates'],
                                                          ascending=(False, True, True, True))
        self.found_decoys = self.found_decoys.reset_index()

        targets_ids = set(
            self.found_decoys['target_id'].to_dict().values())  # get ids of ligands for which decoys were found
        ligands_ids = set(
            self.found_decoys['ligand_id'].to_dict().values())  # get ids of ligands for which decoys were found
        ligands_decoys_saturation = {target: {ligand_id: 0 for ligand_id in ligands_ids} for target in targets_ids}
        ligands_decoys_descriptors = {target: {ligand_id: [] for ligand_id in ligands_ids} for target in targets_ids}
        saturated = set()
        # n_jobs = n_jobs if n_jobs is not None else os.cpu_count()
        for target in targets_ids:
            decoys_dict = {}
            dest = f'{self.storage}/chosen_decoys_for_{target}.smi'
            with open(dest, 'w') as handle:
                for ix, row in self.found_decoys.iterrows():
                    if row['target_id'] == target:
                        ligand_id = row['ligand_id']
                        if ligand_id not in saturated:
                            decoys_dict[row['zinc']] = row['zinc_smile']
                            ligands_decoys_saturation[target][ligand_id] += 1
                            if ligands_decoys_saturation[target][ligand_id] >= top:
                                saturated.add(ligand_id)
                                ligands_descriptors = row[['HBD', 'HBA', 'rotates', 'weight',
                                                           'logp', 'Tanimoto', 'murco_scaffold']]
                                # alphanumeric = "".join([character for character in ligands_descriptors['zinc_smile']
                                #                         if character.isalpha()])  # get all atoms in smiles
                                # ligands_descriptors['zinc_smile'] = len(alphanumeric)  # amount of atoms in smiles (without H)
                                ligands_decoys_descriptors[target][ligand_id].append(ligands_descriptors)
                                if len(saturated) == len(ligands_ids):
                                    print(f'All ligands saturated')
                                    break
                chembl_decoys = load_smiles(f'{target}_inactive.smi')
                for item in chembl_decoys.items():
                    handle.write(f'{item[1]}\t{item[0]}\n')
                for item in decoys_dict.items():
                    handle.write(f'{item[1]}\t{item[0]}\n')
            with open(os.path.join(self.storage, f"saturation.pck"), 'wb') as handle:
                pickle.dump(ligands_decoys_saturation, handle)
            with open(os.path.join(self.storage, f"descriptors.pck"), 'wb') as handle:
                pickle.dump(ligands_decoys_descriptors, handle)

            log(log_file, f'Starting calculating distance matrix for {target}...')
            self.parallel_execute(target, dest)
            enrichment.Prepare(target, os.path.join('chembl_smiles', f'{target}_active.smi'), dest)

    def parallel_execute(self, chembl_id, path):
        Execute(chembl_id, path, self.storage, self.log_file)


class Execute:
    def __init__(self, main_name, decoys_path, storage, log_file):
        with open(decoys_path, 'r') as r:
            self.read_file = r.read().split()
            log(log_file, f"Starting {main_name}")
            self.smiles = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smile)) for smile in self.read_file[::2]]
            log(log_file, f"Made fingerprints for {main_name}")
            self.ids = self.read_file[1::2]
            self.output = pd.DataFrame(columns=self.ids)
            for row in range(len(self.ids)):
                self.output.append(pd.Series(name=self.ids[row]))
            log(log_file, f"CSV output template ready for {main_name}")
            Parallel(n_jobs=1, prefer='threads')(delayed(self.parallel_comparison)(row)
                                                 for row in range(len(self.ids)))
            log(log_file, f"Finished for file {storage}/{main_name}_decoy_bias.csv")
            self.output.to_csv(f'{storage}/{main_name}_decoy_bias.csv')
        print(f'Bailing out from {main_name}')

    def parallel_comparison(self, row_index):
        x = self.smiles[row_index]
        for col in range(row_index, len(self.ids)):
            y = self.smiles[col]
            tc = DataStructs.TanimotoSimilarity(x, y)
            self.output.at[self.ids[row_index], self.ids[col]] = round(tc, 3)
            self.output.at[self.ids[col], self.ids[row_index]] = round(tc, 3)


if __name__ == '__main__':
    Prepare(top=50)
