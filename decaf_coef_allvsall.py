import pandas as pd
import os
import argparse

from rdkit import Chem
from joblib import Parallel, delayed
from decaf.toolkits.rd import phar_from_mol
from decaf.utils import similarity


def log(file_path='', string=''):
    if len(file_path) > 0:
        os.popen(f"echo '{string}' >> {file_path} &")
    else:
        print(string)


class Prepare:
    def __init__(self, path='chembl_smiles', suffix='all.smi', storage='TC_matrices', n_jobs=None, log_file=''):
        self.log_file = log_file
        self.files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path, i)) and suffix in i]
        self.storage = storage
        self.smiles_path = path
        chosen_ones = ['CHEMBL228', 'CHEMBL218', 'CHEMBL2842', 'CHEMBL5071', 'CHEMBL225', 'CHEMBL249',
                         'CHEMBL1978','CHEMBL5102', 'CHEMBL211', 'CHEMBL216', 'CHEMBL2007625', 'CHEMBL231',
                         'CHEMBL2998', 'CHEMBL3572', 'CHEMBL1898', 'CHEMBL4302', 'CHEMBL1821', 'CHEMBL2252',
                         'CHEMBL3772', 'CHEMBL3649', 'CHEMBL2002', 'CHEMBL3836', 'CHEMBL4662', 'CHEMBL5485',
                         'CHEMBL2468', 'CHEMBL3018', 'CHEMBL3898', 'CHEMBL1795186', 'CHEMBL3959', 'CHEMBL4835',
                         'CHEMBL3067', 'CHEMBL2693', 'CHEMBL1781', 'CHEMBL1993', 'CHEMBL1907', 'CHEMBL3243',
                         'CHEMBL1287623', 'CHEMBL3623', 'CHEMBL1926488', 'CHEMBL4769', 'CHEMBL3108640',
                         'CHEMBL2176774', 'CHEMBL3988585', 'CHEMBL3784', 'CHEMBL4941', 'CHEMBL4208']
        self.files_filtered = []
        for i in chosen_ones:
            for j in self.files:
                if i in j:
                    self.files_filtered.append(j)
                    break
        log(log_file, f'Number of chosen CHEMBLS is {len(self.files_filtered)} ')
        self.read_files = []
        if not os.path.exists(storage):
            os.makedirs(storage)
        if len(self.log_file) > 0:
            with open(self.log_file, 'w'): pass
        self.load_files()
        self.parallel_execute()

    def load_files(self):
        for f in self.files_filtered:
            main_name = f.split('_')[0]
            if os.path.isfile(os.path.join(self.storage, '{}_ALLvsALL.csv'.format(main_name))):
                log(self.log_file, "{} exists".format(os.path.join(self.storage, '{}_ALLvsALL.csv'.format(main_name))))
            else:
                with open(os.path.join(self.smiles_path, f), 'r') as f_:
                    f_ = f_.read().split()
                    self.read_files.append((f_, main_name))

    def parallel_execute(self, n_jobs=30):
        Parallel(n_jobs=n_jobs)(delayed(Execute)(i[1], i[0], self.storage, self.log_file) for i in self.read_files)


class Execute:
    def __init__(self, main_name, read_file, storage, log_file):
        log(log_file, "Starting {}".format(main_name))
        smiles = [phar_from_mol(Chem.MolFromSmiles(smile)) for smile in read_file[::2]]
        log(log_file, f"Made pharmacophors for {main_name}")
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
                sim_score, sim_cost = similarity(x, y)
                output.loc[ids[row], ids[col]] = round(sim_score, 3)
                output.loc[ids[col], ids[row]] = round(sim_score, 3)
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
parser.add_argument('--storage', default='Decaf_matrices', help="Folder to which save the outputs")
parser.add_argument('--log_file', default='logs.txt', help="Log to file?")

if __name__ == '__main__':
    args = parser.parse_args()
    Prepare(path=args.path, suffix=args.suffix, storage=args.storage, log_file=args.log_file)
