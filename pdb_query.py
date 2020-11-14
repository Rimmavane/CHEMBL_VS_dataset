from utils import log
from urllib.error import HTTPError
import requests
from json import JSONDecodeError


def get_pdbs_from_unicode(uni_id):
    url = 'https://search.rcsb.org/rcsbsearch/v1/query?'
    queryText = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": f"{uni_id}",
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "UniProt",
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
                    }
                }
            ]
        },
        "return_type": "polymer_entity"
    }
    resp = requests.post(url, json=queryText)
    if resp.ok:
        try:
            result_set = resp.json()['result_set']
            results = [i['identifier'].split('_')[0] for i in result_set]
            resp.close()
            return results
        except JSONDecodeError:
            log(f'Error parsing PDBs for {uni_id}, error code {resp.status_code}')
            return []
    else:
        log(f"Failed to download PDBs for UNIPROT ID {uni_id}")
        return []


def get_pdb_fasta(pdb_id):
    counter = 0
    while True:
        if counter == 5:
            log(f'Stopped trying to download smiles for {pdb_id} after {counter} tries.')
            return ''
        try:
            counter += 1
            url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1'
            resp = requests.get(url)
            if resp.ok:
                fasta = resp.json()['entity_poly']['pdbx_seq_one_letter_code_can']
                resp.close()
                return fasta
            else:
                log(f'Failed to fetch PDB sequence for {pdb_id} with status code {resp.status_code}')
                return ''
        except HTTPError:
            log(f"Could not download ligands PDB sequence for {pdb_id}")
            return ''
        except KeyError:
            log(f"No structure found in PDB {pdb_id}")
            return ''


def check_if_pdb_xray(pdb_code):
    queryText = {
        "query": "query($id: String!){entry(entry_id:$id){exptl{method}}}",
        "variables": {"id": f"{pdb_code}"}
    }
    url = "https://data.rcsb.org/graphql"
    resp = requests.post(url, json=queryText)
    data = resp.json()['data']['entry']
    try:
        data = data['exptl'][0]['method']
        resp.close()
        if 'x-ray' in data.lower() or 'xray' in data.lower():
            return True
    except TypeError:
        return False
    return False


def get_pdb_structure(pdb_id):
    full_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(full_url)
    data = resp.text
    resp.close()
    return data

def get_ligands_for_PDB(pdb_code):
    results = dict()
    try:
        url = "https://data.rcsb.org/graphql"
        queryText = dict(query="""query ($id: String!) {
                                      entry(entry_id:$id) {
                                        nonpolymer_entities {
                                          nonpolymer_comp {
                                            chem_comp {
                                              id
                                            }
                                            pdbx_chem_comp_descriptor{
                                              descriptor
                                              type
                                              program
                                            }
                                          }
                                        }
                                      }
                                    }"""
                         , variables={"id": f"{pdb_code}"})
        resp = requests.post(url, json=queryText)
        assert resp.json()['data']['entry']['nonpolymer_entities']
        ligands = resp.json()['data']['entry']['nonpolymer_entities']
        for i in ligands:
            lig_id = i['nonpolymer_comp']['chem_comp']['id']
            lig_smi = [p['descriptor'] for p in i['nonpolymer_comp']['pdbx_chem_comp_descriptor'] if
                       (p['type'] == 'SMILES_CANONICAL' and p['program'] == 'OpenEye OEToolkits')][0]
            results[lig_id] = lig_smi
        resp.close()
        return results
    except AssertionError:
        log(f'No ligands or query problem for ligands for {pdb_code}')
        return results
