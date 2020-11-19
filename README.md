# CHEMBL_VS_dataset
Code in this repository allows to create a data set basing on data present in ChEMBL. 
Resulting data set contains carefully selected protein targets that fulfilled several conditions, 
including both theoretical and physico-chemcial features.  Each target also has assigned active/inactive 
sets and tools for searching new decoys for actives are provided.

## Requirements
- BLAST+ package for BLAST filtering step.
- Packages listed in requirements.txt


## What pipeline contains

The pipeline consists of few main steps:
1. Initial filter - here targets are being filtered basing on features such as having 
IC50/Kd/Ki activities, existence of PDB structure (X-ray) for each target, checking if ChEMBL targets 
UniProt IDs are present in DUD-E or DEKOIS databases.
Actives and Inactives sets are also created here from all activities for given target, basing on activity threshold.
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

After filtered Actives sets are created, a decoy search for them can be conducted. 
Decoy searching script is separate from the pipeline so it could be used with 
job scheduling programs for spreading computations among more CPUs.

Found decoys than can be filtered to choose N best ones 
(basing on decoy Tanimoto similarity towards ligand) for each ligand.

## Usage

Sample data is available for download here: https://drive.google.com/file/d/1FuJmz8I2rr7yjaOW6ZlOiv42FhiCSILN/view?usp=sharing

Other data can be used as well but it should contain same named columns as in sample data.

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
- --chembl_source   - Accpts single ChEMBL target id (ex. CHEMBL3784) or path to CSV file with 'ChEMBL ID' column form which IDs are retrieved - 
    ligands smi file is found basing on this id, the script searches for <chembl_id>_filtered_active.smi file (or files)
- --chembl_smiles   - folder in which to search for <chembl_id>_filtered_active.smi file, 
    if it don't exists it will be created. By default same as CHEMBL_SMILES_FOLDER parameter in *path_and_settings.py* script
- --decoys_path     - folder with .smi files or path to a single .smi file - if folder given, all .smi files found in 
    this folder will be considered as files with potential decoys
- --output_folder   - where to save results of decoy searching, there will be CSV file 
    for each file found in --decoys_path
- --max_decoys      - how many potential decoys can be assigned to a single ligand in each potential decoys file,
    if there is less files in which the decoys are searched for
- --ignore_pickled  - if found pre-calculated descriptor values for ChEMBL targets ligands 
    in output folder, calculate them again anyway 
--randomize_ligands - randomizing ligands order helps in more even distribution of decoys that satisfies threshold 
    for many ligands

Descriptors thresholds and LogP adaptive thresholding for decoy searching can be set in *paths_and_settings.py* script.

```bash
python3 decoy_searching.py --chembl_source TARGET_ID --chembl_smiles CHEMBL_SMILES_FOLDER --decoys_path FOLDER_WITH_SMI_FILES --output_folder OUTPUT_FOLDER --max_decoys MAX_DECOYS_PER_LIGAND_PER_FILE --ignore_pickled True/False --randomize_ligands True/False
```
Found decoys are saved to files with 'differences.csv' suffix.

Note: Clear or change output_folder if you are going to search decoys for same ChEMBL targets again (ex. with different thresholds),
the results are concatenated to existing CSV files, so if certain target has results from previous searches, than it might cause that
a single decoy would be assigned to more than one ligand

Found decoys should be post-processed with *get_best_found_decoys.py* script, which can be run both from console or by *run_pipeline.py*script.
Parameters accepted during 
- --found_decoys_folder - should be a path to folder with files with 'differences.csv' suffix. 
- --found_decoys_csv    - if not given, is being set to FOUND_DECOYS_CSV/all_found_decoys.csv (from *paths_and_settings*), 
    if a directory path is given, it is used instead of FOUND_DECOYS_CSV, if full path is given it is used as is.
- --decoys_per_ligand   - number of how many best decoys should be chosen for a ligand from all decoys found for that ligand

From console it can be run with command:
```bash
python3 get_best_found_decoys.py --found_decoys_folder FOUND_DECOYS_FOLDER --found_decoys_csv OUTPUT_PATH --decoys_per_ligand BEST_DECOYS_PER_LIGAND
```

Decoys .smi files for all targets for which any decoys were found will be saved to CHEMBL_SMILES_FOLDER (set in *paths_and_settings*) as ID_decoys.smi.
Also in folder given by --found_decoys_folder there will be ID_decoys_amount_found_for_ligands.txt files created with info how many decoys were assigned to each ligand for target of that ID.
