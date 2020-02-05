import pandas as pd
import os
from Bio import SeqIO, Align

# from plotting import similarity_plotting
from pdb_query import get_pdbs_from_unicode, get_pdb_fasta, check_if_pdb_xray
import plotly.graph_objects as go


def choose_primary_pdb_for_chembl(csv_path):
    main_table = pd.read_csv(csv_path, index_col=0)
    main_table['main_PDB_structure'] = ''
    for index, row in main_table.iterrows():
        pdbs = row['PDB_entry'].split()
        fasta = [get_pdb_fasta(i) for i in pdbs]
        best = ''
        for seq in fasta:
            if len(seq) > len(best):
                best = seq
        if best != '':
            best_pdb = pdbs[fasta.index(best)]
            print(best_pdb)
            main_table.at[index, 'main_PDB_structure'] = best_pdb
    main_table.to_csv(csv_path)


def fetch_fastas_for_chembl(csv_path):
    main_table = pd.read_csv(csv_path, index_col=0)
    main_table = main_table[main_table['In_PDBBIND'] == 0]
    print(len(main_table))
    with open(f'{csv_path.strip(".csv").split("/")[-1]}_fastas.txt', 'w') as handle:
        for index, row in main_table.iterrows():
            chembl_name = row['ChEMBL ID']
            pdb_name = row['main_PDB_structure']
            fasta = get_pdb_fasta(row['main_PDB_structure'])
            handle.write(f'>{chembl_name}-{pdb_name}\n{fasta}\n')


def fetch_fastas_for_pdbbind(pdbbind_path):
    pdb_ids = read_pdbbind_ids(pdbbind_path)
    with open('fastas_from_pdbbind.txt', 'w') as handle:
        pdbs_dict = dict()
        for pdb in pdb_ids:
            try:
                if check_if_pdb_xray(pdb):
                    pdb_fasta = get_pdb_fasta(pdb)
                    pdbs_dict[pdb] = pdb_fasta.strip('"')
                    print((pdb_ids.index(pdb) + 1) / len(pdb))
            except:
                print('Not an xray: ', pdb)
        pdbs_dict2 = dict()
        for key, value in pdbs_dict.items():
            if value not in list(pdbs_dict2.values()):
                pdbs_dict2[key] = value
        for key, value in pdbs_dict2.items():
            handle.write(f'>{key}\n{value}\n')


def fetch_fastas_for_DEKOIS(dekois_dir_path):
    dekois_ligands, _ = DEKOIS_uniID_from_folder(dekois_dir_path)
    dekois_ligands_uni_ids = set()
    for i in dekois_ligands:
        try:
            for k in dekois_ligands[i]['Uniprot_ID'].split():
                dekois_ligands_uni_ids.add(k)
        except KeyError:
            print(f'No Uniprot ID for {i}')
    with open('fastas_from_dekois.txt', 'w') as handle:
        for uniprot_id in dekois_ligands_uni_ids:
            pdb_ids = get_pdbs_from_unicode(uniprot_id)
            pdbs_dict = dict()
            for pdb in pdb_ids:
                try:
                    if check_if_pdb_xray(pdb):
                        pdb_fasta = get_pdb_fasta(pdb)
                        pdbs_dict[uniprot_id + '-' + pdb] = pdb_fasta.strip('"')
                        print((pdb_ids.index(pdb) + 1) / len(pdb))
                except:
                    print('Not an xray: ', pdb)
            pdbs_dict2 = dict()
            for key, value in pdbs_dict.items():
                if value not in list(pdbs_dict2.values()):
                    print(key)
                    pdbs_dict2[key] = value
            for key, value in pdbs_dict2.items():
                handle.write(f'>{key}\n{value}\n')


def fetch_fastas_for_DUDE(dude_dir_path):
    dude_uni_ids = DUDE_uniID_from_folder(dude_dir_path)
    with open('fastas_from_dude.txt', 'w') as handle:
        for uniprot_id in dude_uni_ids:
            pdb_ids = get_pdbs_from_unicode(uniprot_id)
            pdbs_dict = dict()
            for pdb in pdb_ids:
                try:
                    if check_if_pdb_xray(pdb):
                        pdb_fasta = get_pdb_fasta(pdb)
                        pdbs_dict[uniprot_id + '-' + pdb] = pdb_fasta.strip('"')
                        print((pdb_ids.index(pdb) + 1) / len(pdb))
                except:
                    print('Not an xray: ', pdb)
            pdbs_dict2 = dict()
            for key, value in pdbs_dict.items():
                if value not in list(pdbs_dict2.values()):
                    print(key)
                    pdbs_dict2[key] = value
            for key, value in pdbs_dict2.items():
                handle.write(f'>{key}\n{value}\n')


#
# def get_targets_fastas(uniprot_id_list, save_path):
#     BASE = 'http://www.uniprot.org'
#     KB_ENDPOINT = '/uniprot/'
#
#     fastas = set()
#     for uniprot_id in uniprot_id_list:
#         payload = {'query': f'id:"{uniprot_id}"',
#                    'format': 'fasta'}
#         result = requests.get(BASE + KB_ENDPOINT, params=payload)
#         if result.ok:
#             fastas.add(result.text)
#         else:
#             print('Something went wrong ', result.status_code)
#     with open(save_path, 'w') as handle:
#         for fasta in fastas:
#             handle.write(fasta)
#     return fastas


def DUDE_uniID_from_folder(dude_folder):
    dude_uni_ids = set()
    for dirpath, dirnames, filenames in os.walk(dude_folder):
        for item in filenames:  # loop through items in dir
            if item.endswith('uniprot.txt'):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of file
                with open(file_name, 'r') as uni_file:
                    for line in uni_file:
                        dude_uni_ids.add(line.strip())
    return dude_uni_ids


def DEKOIS_uniID_from_folder(dir_name, extension='.sdf',
                             search_terms=['Name', 'Uniprot_ID', 'IC50_in_nM', 'Ki_in_nM', 'Kd_in_nM']):
    decoy_dir = dir_name + '/decoys'
    ligand_dir = dir_name + '/ligands'
    uniprot_ligands_ids = dict()
    uniprot_decoys_ids = dict()
    for dirpath, dirnames, filenames in os.walk(ligand_dir):
        for item in filenames:  # loop through items in dir
            if item.endswith(extension):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of files
                with open(file_name, 'r') as sdf_file:
                    for line in sdf_file:
                        for search_term in search_terms:
                            if search_term in line:
                                if search_term == 'Name':
                                    next_line = sdf_file.readline().strip()
                                    uniprot_ligands_ids[next_line] = dict()
                                else:
                                    uniprot_ligands_ids[next_line][search_term] = sdf_file.readline().strip()
    for dirpath, dirnames, filenames in os.walk(decoy_dir):
        for item in filenames:  # loop through items in dir
            if item.endswith(extension):  # check for ".sdf" extension
                file_name = os.path.abspath(os.path.join(dirpath, item))  # get full path of files
                with open(file_name, 'r') as sdf_file:
                    for line in sdf_file:
                        for search_term in search_terms:
                            if search_term in line:
                                if search_term == 'Name':
                                    next_line = sdf_file.readline().strip()
                                    uniprot_decoys_ids[next_line] = dict()
                                else:
                                    uniprot_decoys_ids[next_line][search_term] = sdf_file.readline().strip()

    with open('dekois_ligands_uniprot.txt', 'w') as handle:
        for target in uniprot_ligands_ids:
            handle.write(target + '\n')
            for item in uniprot_ligands_ids[target]:
                handle.write(item + '\t' + uniprot_ligands_ids[target][item] + '\n')

    with open('dekois_decoys_uniprot.txt', 'w') as handle:
        for target in uniprot_decoys_ids:
            handle.write(target + '\n')
            for item in uniprot_decoys_ids[target]:
                handle.write(item + '\t' + uniprot_ligands_ids[target][item] + '\n')

    return uniprot_ligands_ids, uniprot_decoys_ids


def read_pdbbind_ids(file):
    pdb_ids = []
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '#':
                pdb_ids.append(line.split(' ')[0].upper())
    return pdb_ids


def load_fasta_sequences_from_list(fasta_files_list):
    fasta_sequences = {filename.strip('.txt'): [[k.id, k.seq] for k in SeqIO.parse(open(filename), 'fasta')] for
                       filename in fasta_files_list}
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    # aligner.substitution_matrix = MatrixInfo.blosum62
    fasta_sources = list(fasta_sequences.keys())
    for i in range(len(fasta_sources)):
        for j in range(i, len(fasta_sources)):
            if i != j:
                output_name = '-'.join([fasta_sources[i], fasta_sources[j]])
                cols = [v[0] for v in fasta_sequences[fasta_sources[i]]]
                rows = [v[0] for v in fasta_sequences[fasta_sources[j]]]
                values = pd.DataFrame(0.0, index=rows, columns=cols)
                for seq1 in fasta_sequences[fasta_sources[i]]:
                    for seq2 in fasta_sequences[fasta_sources[j]]:
                        alignments = aligner.align(seq1[1], seq2[1])
                        x = str(alignments[0])
                        score = round(x.count('|') / len(seq1[1]), 4)
                        values.at[seq2[0], seq1[0]] = score
                print(values)
                values.to_csv(output_name + '.csv')


def load_fasta_sequences(fasta_file1, fasta_file2, fasta1_to_fasta2_blast):
    fasta_sequences = {''.join(filename.split('.')[:-1]): [[k.id, k.seq] for k in SeqIO.parse(open(filename), 'fasta')]
                       for filename in [fasta_file1, fasta_file2]}
    blast_results = load_blast_results(fasta1_to_fasta2_blast)
    fasta_sources = list(fasta_sequences.keys())
    output_name = '-'.join([fasta_sources[0], fasta_sources[1]])
    cols = [v[0] for v in fasta_sequences[fasta_sources[0]]]
    rows = [v[0] for v in fasta_sequences[fasta_sources[1]]]
    values = pd.DataFrame(0.0, index=rows, columns=cols)
    for val in blast_results:
        values.at[val[1], val[0]] = val[2]
    values.to_csv(output_name + '.csv')


def load_blast_results(blast_results, usecols=[0, 1, 2]):
    df = pd.read_csv(blast_results, header=None, usecols=usecols)
    values = df.values.tolist()
    return values

'''
nazwa nowego targetu, kod struktury, liczba aktywnych, liczba nieaktywnych ligandów, % najwyższego podobieństwa do targetu DEKOIS, e-value, nazwa targetu DEKOIS, % najwyższego podobieństwa do targetu DUDE, e-value, nazwa targetu DUDE, 

Dla kazdego nowego targetu w df jest tylko jedna linia z najwyzszym podobienstwem z obu baz.
Najlepiej jakby przesłał Pan wyniki jak najszybciej (dzis). To zadanie jest juz mocno przeterminowane...
'''


def make_blast_csv(master_path='actives_number_sampling/targets_after_fingerprint_similarity0_tc0.95.csv', blasts=None):
    if blasts is None:
        blasts = ['blast_similarities/chembl-dekois_blast.txt',
                  'blast_similarities/chembl-dude_blast.txt']
    #load master table
    main_table = pd.read_csv(master_path, index_col=0)
    main_table = main_table[main_table['In_PDBBIND'] == 0].copy()
    output_df = main_table[['ChEMBL ID', 'main_PDB_structure', 'Active_compounds', 'Inactive_compounds']]
    output_file = 'chembl_blast_results.csv'
    for db in blasts:
        #chembl_Id, hit_ID, identity%, evalue
        blast_results = load_blast_results(db, [0, 1, 2, 6, 7, 10])
        db_name = db.split('/')[-1].split('_')[0].split('-')[1]
        output_df[f'identity%_{db_name}'] = ''
        output_df[f'evalue_{db_name}'] = ''
        output_df[f'target_name_{db_name}'] = ''
        output_df[f'query_alignment_length_{db_name}'] = ''
        output_df[f'total_query_length_{db_name}'] = ''
        output_df[f'alignment_to_total_ratio_{db_name}'] = ''
        name = ''
        for chembl_id, db_target, identity, q_start, q_end, evalue in blast_results:
            if name not in chembl_id or name == '':
                query_pdb_smile_len = len(get_pdb_fasta(chembl_id.split('-')[-1]))
                query_alignment_length = q_end - q_start
                name = chembl_id.split('-')[0]
                output_df.at[name, f'identity%_{db_name}'] = identity
                output_df.at[name, f'evalue_{db_name}'] = evalue
                output_df.at[name, f'target_name_{db_name}'] = db_target
                output_df.at[name, f'query_alignment_length_{db_name}'] = query_alignment_length
                output_df.at[name, f'total_query_length_{db_name}'] = query_pdb_smile_len
                output_df.at[name, f'alignment_to_total_ratio_{db_name}'] = query_alignment_length/query_pdb_smile_len

    output_df.to_csv(output_file)


def similarity_plotting(csv_path, threshold, title='Similarity heatmap',
                        output_name='similarity_heatmap', size=None, font_size=None):
    main_table = pd.read_csv(csv_path, index_col=0)
    print(main_table.shape)
    main_table = main_table.loc[:, (main_table != 0).any(axis=0)]
    main_table = main_table[(main_table > threshold).any(axis=1)]
    main_table = main_table.loc[:, (main_table > threshold).any(axis=0)]
    print(main_table.shape)

    ####CHECK THE FIGURE
    fig = go.Figure(data=go.Heatmap(
        z=main_table.values.tolist(),
        y=main_table.index.tolist(),
        x=main_table.columns.tolist(),
        colorscale=[[0.0, "blue"],
                    [0.4, "yellow"],
                    [0.6, "orange"],
                    [0.8, "pink"],
                    [1.0, "red"]]))

    fig.update_layout(title=go.layout.Title(text=title, xref="paper", x=0.5))
    if size is not None:
        fig.update_layout(autosize=False, width=size[0], height=size[1])
    fig.update_layout(
        title=go.layout.Title(text=title, xref="paper", x=0.5,
                              font={'size': int(font_size['size'] * 1.3)} if size is not None else None),
        font=font_size)
    fig.write_html(output_name + '.html')
    fig.show()


# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity0_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity1_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity2_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity3_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity4_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity5_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity6_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity7_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity8_tc0.95.csv')
# chembl_fastas = choose_primary_pdb_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity9_tc0.95.csv')


# fetch_fastas_for_DEKOIS('DEKOIS2.0_library')
# fetch_fastas_for_DUDE('DUDE')
# fetch_fastas_for_pdbbind('pdbbind_ids.txt')
# fetch_fastas_for_chembl('actives_number_sampling/targets_after_fingerprint_similarity0_tc0.95.csv')
# load_fasta_sequences('targets_after_fingerprint_similarity10-tc0.95-fastas.txt', 'fastas_from_dude.txt', 'blast_similarities/chembl-dude_blast.txt')
# load_fasta_sequences('targets_after_fingerprint_similarity10-tc0.95-fastas.txt', 'fastas_from_dekois.txt', 'blast_similarities/chembl-dekois_blast.txt')
# load_fasta_sequences('targets_after_fingerprint_similarity10-tc0.95-fastas.txt', 'fastas_from_pdbbind.txt', 'blast_similarities/chembl-pdbbind_blast.txt')

# similarity_plotting(csv_path='targets_after_fingerprint_similarity10-tc095-fastas-fastas_from_dekois.csv',
#                     threshold=95,
#                     title='Sequence similarity heatmap between CHEMBL targets and DEKOIS ligands.',
#                     output_name='heatmap_chembl_dekois')
#                     #,font_size=dict(size=22))
#
# similarity_plotting(csv_path='targets_after_fingerprint_similarity10-tc095-fastas-fastas_from_dude.csv',
#                     threshold=95,
#                     title='Sequence similarity heatmap between CHEMBL targets and DUDE ligands.',
#                     output_name='heatmap_chembl_dude')
#                     #,font_size=dict(size=22))
#
# similarity_plotting(csv_path='targets_after_fingerprint_similarity10-tc095-fastas-fastas_from_pdbbind.csv',
#                     threshold=95,
#                     title='Sequence similarity heatmap between CHEMBL targets and PDBBIND ligands.',
#                     output_name='heatmap_chembl_pdbbind')
#                     #,font_size=dict(size=22))
make_blast_csv()