import pandas as pd
from os import listdir
from os.path import isdir
from utils import log, load_smiles
from rdkit import DataStructs
from utils import make_fingerprint
from paths_and_settings import *
import argparse


def concat_found_decoys_results(found_decoys_folder_path):
    results = []
    for file in listdir(found_decoys_folder_path):
        if file.endswith("differences.csv"):
            x = pd.read_csv(join(found_decoys_folder_path, file), index_col=0)
            results.append(x)
    all_found_decoys = pd.concat(results).drop_duplicates().reset_index(drop=True)
    all_found_decoys.to_csv(join(found_decoys_folder_path, f'all_decoys_found.csv'))
    return all_found_decoys


def add_tanimoto_score(path_found_decoys):
    found_decoys = pd.read_csv(path_found_decoys, index_col=0)
    total = len(found_decoys)
    log(f'Total rows: {total}')
    tanimotos = []
    dec_smiles = dict()
    act_smiles = dict()
    for ix, row in found_decoys.iterrows():
        if ix % 100000 == 0:
            log(f'Added tanimoto scores to {ix} rows')
        try:
            chembl_smile = act_smiles[row['ligand_id']]
        except KeyError:
            chembl_smile = make_fingerprint(row['ligand_smile'], fp='fp2')
            act_smiles[row['ligand_id']] = chembl_smile
        try:
            decoy_smile = dec_smiles[row['zinc']]
        except KeyError:
            decoy_smile = make_fingerprint(row['zinc_smile'], fp='fp2')
            dec_smiles[row['zinc']] = decoy_smile
        tc = DataStructs.FingerprintSimilarity(chembl_smile, decoy_smile)
        tanimotos.append(tc)
    found_decoys = found_decoys.assign(tanimoto=tanimotos)
    found_decoys.to_csv(path_found_decoys)
    return found_decoys


def sort_by_tanimoto(path_to_decoys_csv):
    found_decoys = pd.read_csv(path_to_decoys_csv, index_col=0)
    found_decoys["tanimoto"] = pd.to_numeric(found_decoys["tanimoto"])
    found_decoys = found_decoys.sort_values('tanimoto', ascending=True)
    found_decoys = found_decoys.reset_index(drop=True)
    found_decoys.to_csv(path_to_decoys_csv)
    return found_decoys


def filter_top_hits(path_to_decoys_csv, top=100):
    found_decoys = pd.read_csv(path_to_decoys_csv, index_col=0)
    targets_ids = set(i for i in found_decoys['target_id'].values.tolist())
    for target in targets_ids:
        log(f'Filtering best {top} found decoys for each ligand for {target}')
        target_df = found_decoys[found_decoys['target_id'] == target].drop_duplicates(subset="zinc_smile")
        ligands = list(load_smiles(join(CHEMBL_SMILES_FOLDER, f'{target}_filtered_active.smi')).keys())
        ligands_found_decoys = {ligand: [] for ligand in ligands}
        ligands_decoys_amount = {ligand: 0 for ligand in ligands}

        for ligand in ligands:
            ligand_decoys = target_df[target_df['ligand_id'] == ligand].head(top)
            found_decoys_list = []
            for ix, row in ligand_decoys.iterrows():
                found_decoys_list.append((row['zinc'], row['zinc_smile']))
            ligands_found_decoys[ligand] = found_decoys_list
            ligands_decoys_amount[ligand] = len(found_decoys_list)

        if sum(list(ligands_decoys_amount.values())) > 0:
            log(f'{sum(list(ligands_decoys_amount.values()))} decoys found in total for {target}')
            with open(join(CHEMBL_SMILES_FOLDER, f'{target}_decoys.smi'), 'w') as handle:
                for ligand in ligands_found_decoys:
                    for decoy_tuple in ligands_found_decoys[ligand]:
                        handle.write(f'{decoy_tuple[1]}\t{decoy_tuple[0]}\n')
        else:
            log(f'No decoys found for target {target}')

        with open(join(FOUND_DECOYS_FOLDER, f'{target}_decoys_amount_found_for_ligands.txt'), 'w') as handle:
            for ligand in ligands_decoys_amount:
                handle.write(f'{ligand}\t{ligands_decoys_amount[ligand]}\n')

parser = argparse.ArgumentParser()
parser.add_argument('--found_decoys_folder', help="Folder with files with 'differences.csv' suffixes")
parser.add_argument('--found_decoys_csv', default='', help="Output CSV path")
parser.add_argument('--decoys_per_ligand', default=100, help="How many best decoys should be assigned per ligand")

if __name__ == '__main__':
    args = parser.parse_args()
    FOUND_DECOYS_STORAGE = args.found_decoys_folder
    FOUND_DECOYS_CSV = args.found_decoys_csv
    if FOUND_DECOYS_CSV == '':
        FOUND_DECOYS_CSV = join(FOUND_DECOYS_STORAGE, 'all_decoys_found.csv')
    elif isdir(FOUND_DECOYS_CSV):
        FOUND_DECOYS_CSV = join(FOUND_DECOYS_CSV, 'all_decoys_found.csv')

    _ = concat_found_decoys_results(FOUND_DECOYS_STORAGE)
    _ = add_tanimoto_score(FOUND_DECOYS_CSV)
    _ = sort_by_tanimoto(FOUND_DECOYS_CSV)
    filter_top_hits(FOUND_DECOYS_CSV, top=int(args.decoys_per_ligand))
