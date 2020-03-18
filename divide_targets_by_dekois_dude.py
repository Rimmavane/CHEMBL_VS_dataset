import pandas as pd
import os


def sample_by_active_compounds(list_of_values):
    if not os.path.exists('actives_number_sampling'):
        os.makedirs('actives_number_sampling')
    targets = pd.read_csv('master_table.csv', sep=',', index_col=0)
    print(f'Number of targets at the beginning: {len(targets)}')
    new = targets[pd.notna(targets['DEKOIS_ID']) | pd.notna(targets['DUDE_ID'])]
    print(f'Targets with coverage: {len(new)}')
    new.to_csv(f"actives_number_sampling/targets_with_coverage.csv")
    for val in list_of_values:
        new = targets[pd.isna(targets['DEKOIS_ID']) & pd.isna(targets['DUDE_ID'])]
        new = new.sort_values(by=['Active_compounds'], ascending=False)
        new = new[new['Active_compounds'] >= val]
        new = new.drop(labels = ['DEKOIS_ID', 'DEKOIS_actives',
       'DEKOIS_actives_threshold', 'DEKOIS_decoys', 'DUDE_ID',
       'Active_in_DUDE', 'Active_in_DUDE_threshold', 'Inactive_in_DUDE',
       'Inactive_in_DUDE_threshold'], axis=1)
        print('Targets without coverage:', len(new))
        new.to_csv(f"actives_number_sampling/targets_without_coverage_{val}.csv")

# targets = pd.read_csv('master_table.csv', sep=',', index_col=0)
# print(f'Number of targets at the beginning: {len(targets)}')
# new = targets[pd.notna(targets['DEKOIS_ID']) | pd.notna(targets['DUDE_ID'])]
# print(f'Number of CHEMBL targets found in DEKOIS or DUDE: {len(new)}')
# new.to_csv("targets_with_coverage.csv")
# ########################################
# new = targets[pd.isna(targets['DEKOIS_ID']) & pd.isna(targets['DUDE_ID'])]
# new = new.sort_values(by=['Active_compounds'], ascending=False)
# print(f'Number of CHEMBL targets NOT found in DEKOIS or DUDE: {len(new)}')
# new = new[new['Active_compounds'] >= 50]
# print(f'Number of CHEMBL targets with more than 50 actives: {len(new)}')
# new = new.drop(labels=['DEKOIS_ID', 'DEKOIS_actives',
#                        'DEKOIS_actives_threshold', 'DEKOIS_decoys', 'DUDE_ID',
#                        'Active_in_DUDE', 'Active_in_DUDE_threshold', 'Inactive_in_DUDE',
#                        'Inactive_in_DUDE_threshold'], axis=1)
# new.to_csv("targets_without_coverage.csv")

sample_by_active_compounds([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
