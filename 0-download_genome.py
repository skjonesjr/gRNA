# %% [markdown]
# # Download and unzip human genome

from pathlib import Path
from datetime import datetime
import urllib.request
import os, gzip, shutil

import tqdm


with open('.env') as f:
    CONFIG = dict([line.strip().split('=', maxsplit=1) for line in f.readlines()])
RAW = Path(CONFIG['RAW'])


def download_and_extract_gz(url, output_dir):
    # Download the gzipped file
    file_name = url.split('/')[-1]
    file_path = os.path.join(output_dir, file_name)
    urllib.request.urlretrieve(url, file_path)

    # Extract the contents
    with gzip.open(file_path, 'rb') as f_in:
        with open(file_path[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Remove the gzipped file
    os.remove(file_path)


url = 'https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{}.fa.gz'


# Should have about 3.5 GB of free space
output_dir = RAW / 'human_genome'
output_dir.mkdir(exist_ok=True, parents=True)

with open(output_dir / 'README', 'w') as f:
    f.write(f'Date: {datetime.now().strftime("%Y-%m-%d")}\n')
    f.write(f'Source: {url}\n')

for chrom in tqdm.tqdm(list(range(1, 23)) + list('XY')):
    # Takes a long time to extract
    download_and_extract_gz(url.format(chrom), output_dir)
