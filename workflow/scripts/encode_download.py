import pdb
import os
from pathlib import Path
import wget
import pandas as pd
from src.utils.io import read_data_matrix


df = read_data_matrix(snakemake.input['data_matrix'])
encode_id = snakemake.wildcards['encode_id']

if encode_id in df.index:
    url = 'file_url'
else:
    df = df.set_index('paired_acc')
    url = 'paired_url'

row = df.loc[encode_id]

folder = Path(snakemake.output[0]).parent

filename = str(folder / f'{row.name}.fastq.gz')
wget.download(row[url], out=filename)
os.system(f'gzip -d {filename}')
