from os.path import join
from os import getcwd
from utils import create_folder_if_not_existent
import time

# IF LOGGING SET TO TRUE, COMMUNICATES WILL BE REDIRECTED TO A FILE, OTHERWISE THEY WILL BE PRINTED IN CONSOLE

# MAIN PATH FROM WHICH THE PIPELINE WILL BE WORKING
PROJECT_HOME = getcwd()
MAIN_CSV_NAME = join(PROJECT_HOME, 'master_table.csv')

LOGGING = False
LOG_FILE = ''
if LOGGING:
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y', t)
    LOG_FILE = join(PROJECT_HOME, f'log_{timestamp}.txt')
    with open(LOG_FILE, 'w') as _:
        pass

# REQUIRED INPUT DATA
RAW_DATA_FOLDER = join(PROJECT_HOME, 'raw_data')
RAW_TARGETS = join(RAW_DATA_FOLDER, 'CHEMBL25-targets.csv')                 # CSV from with targets will be taken
RAW_IC50_ACTIVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_ic50.csv')       # CSV with IC50 activities for targets
RAW_KI_ACTIVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_ki.csv')           # CSV with Ki activities for targets
RAW_KD_ACTVITY = join(RAW_DATA_FOLDER, 'CHEMBL25-single_kd.csv')            # CSV with Kd activities for targets
RAW_UNIPROT_3D_IDS = join(RAW_DATA_FOLDER, 'has_3d_pdb_acc.csv')            # file with uniprot ID's for which any X-ray PDB structure exists

DEKOIS_PATH = join(PROJECT_HOME, 'DEKOIS2.0_library')                       # path to folder containing DEKOIS database
DUDE_PATH = join(PROJECT_HOME, 'DUDE')                                      # path to folder containing DUD-E database


# OUTPUTS DURING PIPELINE
STANDARDS_FOLDER = join(PROJECT_HOME, 'standards_csv')                      # processed activities CSV's will be stored here
create_folder_if_not_existent(STANDARDS_FOLDER)

COMPOUND_SAMPLING_FOLDER = join(PROJECT_HOME, 'actives_number_sampling')    # folder in which results of filtering by actives number are stored
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
LIGAND_WEIGHT_LOWER_THRESHOLD = 100  # int
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
USE_JOBLIB = False
JOBLIB_WORKERS = 4


# DECOY SEARCHING
# logp difference is calculated by abs(difference)-(1-(abs(chembl_ligands_logp)/10)), so for ligands with very low logp any hits were still possible
# and it narrowed possible matches for ligands with very high logp
# (because 20% from 0 would give 0, so that would be very rough threshold for logp, on the other hand 20% from 5 is less than 20% from 8 so the bigger the value the less added value it gets)
THRESHOLDS_DICT = {'HBD': 0, 'HBA': 0, 'rotates': 0, 'weight': 0.15, 'logp': 0.2, 'murco_scaffold': 0.7}  # sample values
