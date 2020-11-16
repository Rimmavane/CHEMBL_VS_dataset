from os import listdir
from os.path import isfile
from urllib.error import URLError
import urllib.request
import argparse
from paths_and_settings import *


def download_zinc(file_path, dest_path=join(RAW_DATA_FOLDER, 'ZINC'), overwrite=False):
    create_folder_if_not_existent(dest_path)
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
                    urllib.request.urlretrieve(url, join(dest_path, name))
                print(f'DOWNLOADED: {url}')
            except URLError:
                print(f'FATAL ERROR FOR {url}')


def _str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--zinc_file', help="Folder with links to ZINC resources")
parser.add_argument('-o', '--output', default=f'{join(RAW_DATA_FOLDER, "ZINC")}',
                    help="Folder where downloaded files will be stored")
parser.add_argument('--overwrite', type=_str2bool, default=False,
                    help="If change to True it will be checked if resource is already in output "
                         "folder to avoid downloading it again")

if __name__ == '__main__':
    args = parser.parse_args()
    zinc_file = args.zinc_file
    output = args.output
    overwrite = args.overwrite
    download_zinc(file_path=zinc_file, dest_path=output, overwrite=overwrite)
