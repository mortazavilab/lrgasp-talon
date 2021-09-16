import pandas as pd


def read_data_matrix(path):
    df = pd.read_csv(path, sep='\t')
    df = df[df['file_contents'] == 'reads']
    df['sample'] = df['sample'].str.replace('_', '')
    return df.set_index('file_acc')
