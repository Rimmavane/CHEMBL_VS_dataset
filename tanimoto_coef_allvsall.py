import pandas as pd
import os

from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from joblib import Parallel, delayed
import argparse


def log(file_path='', string=''):
    if len(file_path) > 0:
        os.popen(f"echo '{string}' >> {file_path} &")
    else:
        print(string)


class Prepare:
    def __init__(self, path='chembl_smiles', suffix='active.smi', storage='TC_matrices', n_jobs=None, log_file=''):
        self.log_file = log_file
        self.files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and suffix in i]
        self.storage = storage
        self.smiles_path = path
        self.read_files = []
        n_jobs = n_jobs if n_jobs is not None else os.cpu_count()-1
        if not os.path.exists(storage):
            os.makedirs(storage)
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'): pass
        log(self.log_file, f"Availible cores: {n_jobs}")
        self.load_files()
        self.parallel_execute(n_jobs)

    def load_files(self):
        for f in self.files:
            main_name = f.split('_')[0]
            if os.path.isfile(os.path.join(self.storage, '{}_ALLvsALL.csv'.format(main_name))):
                log(self.log_file, "{} exists".format(os.path.join(self.storage, '{}_ALLvsALL.csv'.format(main_name))))
            else:
                with open(os.path.join(self.smiles_path, f), 'r') as f_:
                    f_ = f_.read().split()
                    self.read_files.append((f_, main_name))

    def parallel_execute(self, n_jobs):
        Parallel(n_jobs=n_jobs)(delayed(Execute)(i[1], i[0], self.storage, self.log_file) for i in self.read_files)


class Execute:
    def __init__(self, main_name, read_file, storage, log_file):
        log(log_file, "Starting {}".format(main_name))
        smiles = [FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smile)) for smile in read_file[::2]]
        log(log_file, f"Made fingerprints for {main_name}")
        ids = read_file[1::2]
        output = pd.DataFrame(columns=ids)
        s = 0
        for row in range(len(ids)):
            output.append(pd.Series(name=ids[row]))
        log(log_file, f"CSV output template ready for {main_name}")
        for row in range(len(ids)):
            x = smiles[row]
            for col in range(row, len(ids)):
                y = smiles[col]
                tc = DataStructs.FingerprintSimilarity(x, y)
                output.loc[ids[row], ids[col]] = round(tc, 3)
                output.loc[ids[col], ids[row]] = round(tc, 3)
            if s % 10 == 0:
                s += 1
                log(log_file, f"{s}/{len(ids)} smiles checked for {main_name}")
                s -= 1
            s += 1
        output.to_csv('{}/{}_ALLvsALL.csv'.format(storage, main_name))
        log(log_file, f"Finished for file {storage}/{main_name}_ALLvsALL.csv")


parser = argparse.ArgumentParser()
parser.add_argument('--path', default='chembl_smiles', help="Folder with CHEMBL smiles files")
parser.add_argument('--suffix', default='all.smi')
parser.add_argument('--storage', default='TC_matrices', help="Folder to which save the outputs")
parser.add_argument('--log_file', default='logs.txt', help="Log to file?")

if __name__ == '__main__':
    args = parser.parse_args()
    Prepare(path=args.path, suffix=args.suffix, storage=args.storage, log_file=args.log_file)
