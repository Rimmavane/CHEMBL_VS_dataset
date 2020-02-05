import os
import pandas as pd
from progress.bar import Bar
from tanimoto_coef_allvsall import log
from joblib import Parallel, delayed
from smilite import get_zinc_smile
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


def create_folder_if_not_existent(path):
    if not os.path.exists(path):
        os.makedirs(path)


def load_smiles(filepath):
    smiles_dict = dict()
    # if path to file provided
    try:
        with open(filepath, 'r') as smiles_file:
            for line in smiles_file:
                line = line.split()
                smiles_dict[line[1]] = line[0]
    # if only chembl name provided
    except FileNotFoundError:
        file = \
        [i for i in os.listdir('chembl_smiles') if os.path.isfile(os.path.join('chembl_smiles', i)) and filepath in i][0]
        with open(f'chembl_smiles/{file}', 'r') as smiles_file:
            for line in smiles_file:
                line = line.split()
                smiles_dict[line[1]] = line[0]
    return smiles_dict


class Prepare:
    def __init__(self, chosen_chembl, path_found_decoys='found_decoys', top=65,
                 storage='decoys_matrices', n_jobs=None, log_file=''):
        self.log_file = log_file
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'): pass

        log(log_file, 'Beginning loading found decoys...')
        self.found_decoys = self.load_decoys(path_found_decoys, how_many=4)

        log(log_file, 'Beginning filtering found decoys...')
        self.all_found_decoys = {i: self.filter_decoys_by_chembl_id(i) for i in chosen_chembl}
        log(log_file, 'Beginning loading CHEMBL decoys...')
        self.chembl_decoys = {i: load_smiles(f'{i}_inactive.smi') for i in chosen_chembl}

        self.storage = storage
        create_folder_if_not_existent(storage)

        # self.read_files = []
        n_jobs = n_jobs if n_jobs is not None else os.cpu_count()
        log(log_file, 'Starting calculating distance matrix...')
        for chembl_id in chosen_chembl:
            p = self.get_top_decoys(chembl_id, top)
            self.parallel_execute(n_jobs, chembl_id, p)

    def load_decoys(self, found_decoys_folder, how_many):
        gatherer = dict()
        names = set([i.split('.smi')[0] for i in os.listdir(found_decoys_folder) if
                     os.path.isfile(os.path.join(found_decoys_folder, i)) and 'values' not in i])
        search_patterns = [f'_{i}_' for i in range(how_many)]
        self.threshold = [i for i in range(how_many)]
        for file_name in names:
            file_csv = dict()
            for pattern in search_patterns:
                ix = search_patterns.index(pattern)
                file_threshold_csv = None
                file = [i for i in os.listdir(found_decoys_folder)
                        if os.path.isfile(os.path.join(found_decoys_folder, i)) and pattern in i
                        and file_name in i and 'values' not in i]
                file = file[0] if len(file) > 0 else None
                if file is not None:
                    file_threshold_csv = pd.read_csv(os.path.join(found_decoys_folder, file), index_col=0)
                    file_threshold_csv['file'] = file
                    #     new_csv = pd.read_csv(os.path.join(found_decoys_folder, file), index_col=0)
                    #     new_csv['file'] = file
                    #     file_threshold_csv = pd.concat([file_threshold_csv, new_csv]).drop_duplicates().reset_index(drop=True)

                file_csv[ix] = file_threshold_csv
            gatherer[file_name] = file_csv
        return gatherer

    def filter_decoys_by_chembl_id(self, chembl_id):
        filtered_decoys = dict()
        # file names
        for file in self.found_decoys:
            filtered_decoys[file] = dict()
            # thresholds
            for threshold in self.found_decoys[file]:
                if self.found_decoys[file][threshold] is None:
                    filtered_decoys[file][threshold] = None
                else:
                    filtered_decoys[file][threshold] = self.found_decoys[file][threshold][
                        self.found_decoys[file][threshold]['target_id'] == chembl_id]
                    filtered_decoys[file][threshold] = filtered_decoys[file][threshold].sort_values(
                        by=['HBD', 'HBA', 'rotates', 'weight', 'logp'], axis=0)
        return filtered_decoys

    def get_top_decoys(self, chembl_id, top=200):
        decoys = self.chembl_decoys[chembl_id].copy()
        smiles = set(decoys.values())
        found_decoys_chembl_id = self.all_found_decoys[chembl_id]
        if len(decoys) < top:
            for i in self.threshold:
                log(self.log_file, "Current smiles formatted {} before threshold {}".format(len(decoys), i))
                temp = None
                for file in found_decoys_chembl_id:
                    decoys_from_file_with_th = found_decoys_chembl_id[file][i]
                    if decoys_from_file_with_th is None:
                        pass
                    else:
                        temp = decoys_from_file_with_th if temp is None else pd.concat([temp, decoys_from_file_with_th])
                if temp is None:
                    pass
                else:
                    temp = temp.sample(frac=1).reset_index(drop=True)
                    unique = temp['zinc'].unique().tolist()
                    bar = Bar('Checking ZINCs:', max=len(unique))
                    for zinc in unique:
                        zinc_smi = get_zinc_smile(zinc, backend='zinc15')
                        if zinc_smi not in smiles:
                            decoys[zinc] = zinc_smi
                            smiles.add(zinc_smi)
                            if len(decoys) >= top:
                                log(self.log_file, ' Max number of decoys reached, bailing off early...')
                                break
                        bar.next()
                    bar.finish()
                if len(decoys) >= top:
                    break
        with open(f'{self.storage}/{chembl_id}_filtered_decoys.smi', 'w') as handle:
            for item in decoys.items():
                handle.write(f'{item[1]}\t{item[0]}\n')
        return f'{self.storage}/{chembl_id}_filtered_decoys.smi'

    def parallel_execute(self, n_jobs, chembl_id, path):
        Execute(chembl_id, path, self.storage, self.log_file)


class Execute:
    def __init__(self, main_name, p, storage, log_file):
        with open(p, 'r') as r:
            self.read_file = r.read().split()
            log(log_file, f"Starting {main_name}")
            self.smiles = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smile)) for smile in self.read_file[::2]]
            log(log_file, f"Made fingerprints for {main_name}")
            self.ids = self.read_file[1::2]
            self.output = pd.DataFrame(columns=self.ids)
            for row in range(len(self.ids)):
                self.output.append(pd.Series(name=self.ids[row]))
            log(log_file, f"CSV output template ready for {main_name}")
            Parallel(n_jobs=1, prefer='threads')(delayed(self.parallel_comparison)(row)
                                                 for row in range(len(self.ids)))
            self.output.to_csv('{}/{}_decoy_bias.csv'.format(storage, main_name))
            log(log_file, f"Finished for file {storage}/{main_name}_ALLvsALL.csv")

    def parallel_comparison(self, row_index):
        x = self.smiles[row_index]
        for col in range(row_index, len(self.ids)):
            y = self.smiles[col]
            tc = DataStructs.FingerprintSimilarity(x, y)
            self.output.at[self.ids[row_index], self.ids[col]] = round(tc, 3)
            self.output.at[self.ids[col], self.ids[row_index]] = round(tc, 3)


if __name__ == '__main__':
    Prepare(chosen_chembl=['CHEMBL3623'], top=750)
