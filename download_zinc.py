import wget
from os import listdir
from os.path import isfile, join

def download_zinc(file_path, dest_path, overwrite = False):
    downloaded_files = [f for f in listdir(dest_path) if isfile(join(dest_path, f))]
    with open(file_path, 'r') as handle:
        for line in handle:
            url = line.strip()
            name = url.split('/')[-1]
            exists = 0
            try:
                if overwrite:
                    exists = 0
                else:
                    for file_name in downloaded_files:  # check if already downloaded
                        if name in file_name:
                            exists = 1
                            break
                if exists == 0:
                    wget.download(url, f'{dest_path}/{name}')
                print(f'DOWNLOADED: {url}')
            except:
                print(f'FATAL ERROR FOR {url}')


def main():
    download_zinc('D:/Studia/Magisterka/chembl/ZINC-downloader-2D-smi.uri', 'D:/Studia/Magisterka/chembl/raw_data/ZINC')

main()