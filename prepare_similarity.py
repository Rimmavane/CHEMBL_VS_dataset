import pandas as pd
import os
from rdkit import DataStructs
from utils import make_fingerprint, load_smiles,  log
from paths_and_settings import *


class Prepare_Similarity:
    def __init__(self, output_name, path_query, path_subject, storage=SIMILARITY_MATRICES_FOLDER, do_log=LOGGING):
        self.log_file = False
        self.storage = storage
        self.name = output_name
        if not os.path.exists(storage):
            os.makedirs(storage)
        if do_log is True:
            self.log_file = LOG_FILE
        else:
            self.log_file = ""
        self.queries = load_smiles(path_query)
        self.subjects = load_smiles(path_subject)
        self.execute()

    def execute(self):
        queries = self.queries
        subjects = self.subjects
        log(f"Starting {self.name} queries and subjects")
        queries = {i: make_fingerprint(queries[i], fp='fp2') for i in queries}
        subjects = {i: make_fingerprint(subjects[i], fp='fp2') for i in subjects}

        log(f"Made fingerprints for {self.name} queries and subjects")
        output = pd.DataFrame(index=list(queries.keys()), columns=list(subjects.keys()))
        output_path = os.path.join(self.storage, f'{self.name}.csv')
        s = 0
        log(f"CSV output template ready for {self.name} queries and subjects")
        for row in queries:
            x = queries[row]
            for col in subjects:
                y = subjects[col]
                tc = DataStructs.FingerprintSimilarity(x, y)
                output.at[row, col] = round(tc, 3)
            if s % 25 == 0:
                log(f"{s}/{len(queries)} smiles checked for {self.name}")
            s += 1
        output.to_csv(output_path)
        log(f"Finished for file {output_path}")


def prepare_actives_similarirty_matrices(main_csv_path):
    main_csv = pd.read_csv(main_csv_path, index_col=0)
    ids = main_csv['ChEMBL ID'].tolist()
    for target in ids:
        actives_path = os.path.join(CHEMBL_SMILES_FOLDER, f'{target}_active.smi')
        log(f'Starting counting AvA similarity matrix for {target}')
        Prepare_Similarity(output_name=f'{target}_AvA_similarity',
                           path_query=actives_path,
                           path_subject=actives_path)
