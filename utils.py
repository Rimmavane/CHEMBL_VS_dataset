from rdkit.Chem import AllChem, MACCSkeys
from rdkit import Chem
import os
from itertools import groupby


def make_fingerprint(smile: str, fp='fp2'):
    """
    Make fingerprint for given smile
    :param smile: Compound smile
    :return: fingerprint
    """
    try:
        x = Chem.MolFromSmiles(smile)
        if fp == 'fp2':
            x = Chem.RDKFingerprint(x)
        elif fp == 'morg2':
            x = AllChem.GetMorganFingerprintAsBitVect(x, 2, useFeatures=True)
        else:
            x = MACCSkeys.GenMACCSKeys(x)
        return x
    except:
        log(f'COULDNT MAKE FINGERPRINT FROM {smile}')
        return None


def create_folder_if_not_existent(path):
    if not os.path.exists(path):
        os.makedirs(path)


def load_smiles(filepath):
    smiles_dict = dict()
    with open(filepath, 'r') as smiles_file:
        for line in smiles_file:
            line = line.strip().split()
            smiles_dict[line[1]] = line[0]
    return smiles_dict


def log(string=''):
    from paths_and_settings import LOG_FILE
    if len(LOG_FILE) > 0:
        os.popen(f"echo '{string}' >> {LOG_FILE} &")
    else:
        print(string)


def read_fasta(fasta_path):
    """Generator over the fasta files"""
    fasta_file = open(fasta_path, 'r')
    fasta_iter = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))

    for header in fasta_iter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in fasta_iter.__next__())
        yield (headerStr, seq)


def pandas_col_to_file(pandas_col, output_path):
    pandas_col.to_csv(output_path, index=False, header=False)


def get_files_with_suffix(folder_path, suffix='_filtered_active.smi'):
    """
    Find all files with given suffix in specified folder, without recursion
    :param folder_path: Folder we would like to search in
    :param suffix: Suffix of files we want to find
    :return: List of paths to files
    """
    files = [os.path.join(folder_path, i) for i in os.listdir(folder_path) if
             os.path.isfile(os.path.join(folder_path, i)) and suffix in i]
    return files
