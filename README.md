# CHEMBL_VS_dataset
Code in this repository allows to create a data set basing on data present in ChEMBL. Resulting data set contains carefully selected protein targets that fullfilled several conditions, including both theoretical and physico-chemcial features.  Each target also has assigned active/inactive sets and tools for searching new decoys for actives are provided.

## Requirements
- BLAST+ package for BLAST filtering step.
- Packages listed in requirements.txt


## What pipeline contains

The pipeline consists of few main steps:
1. Initial filter - here targets are being filtered basing on features such as having 
IC50/Kd/Ki activities, existance of PDB stucture (X-ray) for each target, checking if ChEMBL targets 
UniProt IDs are present in DUD-E or DEKOIS databases.
Actives and Inactives sets are also created here from all activities for given target, basing on activity thereshold.
2. Actives number filter - the filters are: number of compounds in Active set, existance of ligands (within PDB structures) 
similar to compounds in Actives set (basing on Tanimoto coefficient). At this stage sequences for PDBs are downloaded and a primary 
PDB structure is being chosen basing on length of sequences.
3. BLAST filter - Targets (sequences of primary PDB structures) are BLASTed against active 
compounds from DUD-E and DEKOIS databases. Than all targets containing hits with BLAST similarities 
over given threshold are filtered out.
4. Creating Tanimoto similarity matrices for Actives sets
5. Filtering Actives basing on similarity threshold - if a pair of actives have similarity over 
threshold, the one with worse (higher) minimal activity is discarded. This step reduces analogue 
bias within Actives set - compounds are not structurally similar beyond given threshold.

After filtered Actives sets are created, a decoy serach for them can be conducted. 
Decoy searching script is separate from the pipeline so it could be used with 
job scheduling programs for spreading computations among more CPUs.

## Usage

Scripts *paths_and_settings.py* and *run_pipeline.py* are used as the main part of managing the filtering process.
*paths_and_settings.py* contains all major paths of input data as well as output data. 
It is also where the filtering and decoy searching thresholds can be set by user.

Decoy searching code can be used with any set of smiles files but this package provides code for downloading ZINC files.
To download ZINC files first head to ZINC database website and download file with links to resources of certain tranche.
Than use *download_zinc.py* script, it can be used from console with command:

```bash
python3 download_zinc.py -f ZINC_FILE -o OUTPUT_FOLDER
```

*decoy_searching.py* script can be used from console as well, possible arguments are:
- --chembl_id      - ChEMBL target id (ex. CHEMBL3784) - query smiles file is determined based on this id,
the script searches for <chembl_id>_filtered_active.smi file
- --chembl_smiles  - folder in which to search for <chembl_id>_filtered_active.smi file, if it don't exists it will be created
- --decoys_path    - folder with .smi files - all .smi files found in this folder will be considered as files with potential decoys
- --output_folder  - where to save results of decoy searching, there will be CSV file for each file found in --decoys_path
- --max_decoys     - how many potential decoys can be assigned to a single ligand in each potential decoys file
- --ignore_pickled - if found pre-calculated descriptor values for ChEMBL targets ligands in output folder, calculate them again anyway 

Note: Clear or change output_folder if you are going to search decoys for same ChEMBL targets again (ex. with different thresholds),
the results are concatenated to existing CSV files, so if certain target has results from previous searches, than it might cause that
a single decoy would be assigned to more than one ligand

TODO:
- write more detailed README
- add ligand shuffling in decoy searching
- add code for found decoys post-processing
