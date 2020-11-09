from os.path import join
from os import getcwd
from utils import create_folder_if_not_existent
import time

# IF LOGGING SET TO TRUE, COMMUNICATES WILL BE REDIRECTED TO A FILE, OTHERWISE THEY WILL BE PRINTED IN CONSOLE

# MAIN PATH FROM WHICH THE PIPELINE WILL BE WORKING
PROJECT_HOME = getcwd()
MAIN_CSV_NAME = join(PROJECT_HOME, 'master_table.csv')

LOGGING = True
LOG_FILE = ''
if LOGGING:
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y', t)
    LOG_FILE = join(PROJECT_HOME, f'log_{timestamp}.txt')

# REQUIRED INPUT DATA
RAW_DATA_FOLDER = join(PROJECT_HOME, 'raw_data')
RAW_TARGETS = join(RAW_DATA_FOLDER, 'CHEMBL25-targets.csv')
RAW_IC50_ACTIVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_ic50.csv')
RAW_KI_ACTIVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_ki.csv')
RAW_KD_ACTVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_kd.csv')
RAW_UNIPROT_3D_IDS = join(RAW_DATA_FOLDER, 'has_3d_pdb_acc.csv')

DEKOIS_PATH = join(PROJECT_HOME, 'DEKOIS2.0_library')
DUDE_PATH = join(PROJECT_HOME, 'DUDE')


# OUTPUTS DURING PIPELINE
STANDARDS_FOLDER = join(PROJECT_HOME, 'standards_csv')
create_folder_if_not_existent(STANDARDS_FOLDER)

COMPOUND_SAMPLING_FOLDER = join(PROJECT_HOME, 'actives_number_sampling')
create_folder_if_not_existent(COMPOUND_SAMPLING_FOLDER)

BLAST_MAIN_FOLDER = join(PROJECT_HOME, 'blast_similarities')
create_folder_if_not_existent(BLAST_MAIN_FOLDER)

CHEMBL_SMILES_FOLDER = join(PROJECT_HOME, 'chembl_smiles')
create_folder_if_not_existent(CHEMBL_SMILES_FOLDER)

SIMILARITY_MATRICES_FOLDER = join(PROJECT_HOME, 'similarity_matrices')
create_folder_if_not_existent(SIMILARITY_MATRICES_FOLDER)

# FILTERING ARGUMENTS
INITIAL_FILTER = True
ACTIVITY_VALUE_THRESHOLD = 10  # in nM

SAMPLING_FILTER = True
LIGAND_WEIGHT_LOWER_THRESHOLD = 200  # int
LOWEST_TC_SIMILARITY_BETWEEN_LIGANDS_THRESHOLD = (0.95,)    # a tuple
LOWER_LIMIT_OF_LIGANDS = (50,)  # a tuple

# BLAST ARGUMENTS
BLAST = True
CHOSEN_TC_THRESHOLD = 0.95  # float
CHOSEN_LIGAND_LIMIT = 50    # int
E_VALUE_THRESHOLD = 0.00001   # float
MAX_BLAST_SIMILARITY = 30  # sequences % similarity threshold between targets and best hits from DEKOIS/DUD-E

# CRATING SIMILARITY MATRICES
CREATE_ACTIVES_SIMILARITY_MATRICES = True

# REDUCING ANALOGUE BIAS AMONG ACTIVES
FILTER_ACTIVES = True
ACTIVES_TC_SIMILARITY_THRESHOLD = 0.9  # maximum Tanimoto similarity between compounds in active set
