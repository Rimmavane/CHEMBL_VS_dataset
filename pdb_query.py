import urllib.request
import pypdb

def get_pdbs_from_unicode(uni_id):
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = (f"    \n"
                 f"    <orgPdbQuery>\n"
                 f"    \n"
                 f"    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>\n"
                 f"    \n"
                 f"    <description>Simple query for a list of Uniprot Accession IDs: {uni_id}</description>\n"
                 f"    \n"
                 f"    <accessionIdList>{uni_id}</accessionIdList>\n"
                 f"    \n"
                 f"    </orgPdbQuery>\n"
                 f"    ").encode('utf-8')
    req = urllib.request.Request(url)

    with urllib.request.urlopen(req, data=queryText) as f:
        resp = f.read().decode('utf-8').split()
        resp = [i[:4] for i in resp]
        if 'null' in resp or '' in resp:
            return []
        return resp


def get_pdb_fasta(pdb_id):
    while True:
        counter = 0
        if counter == 200:
            raise Exception(f'Stopped trying to download smiles for {pdb_id} after {counter} tries.')
        try:
            with urllib.request.urlopen(f'http://www.rcsb.org/pdb/rest/customReport.csv?pdbids={pdb_id}&customReportColumns=sequence&format=csv&service=wsfile') as f:
                resp = f.read()
                fastas = resp.decode("utf-8").replace('\n', ',').split(',')[5::3]
                best = ''
                for fasta in fastas:
                    if len(fasta) > len(best):
                        best = fasta
                return best.strip('"')
        except:
            counter += 1
            continue


def check_if_pdb_xray(pdb_code):
    if pypdb.get_entity_info(pdb_code)['Method']['@name'] == 'xray':
        return True
    return False
