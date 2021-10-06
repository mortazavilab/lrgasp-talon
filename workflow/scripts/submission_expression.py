import pdb
import pandas as pd

df = pd.read_csv(snakemake.input['abundance'], sep='\t')
samples = df.columns[11:].tolist()
df = df[['annot_transcript_id'] + samples]

renames = {i: i.split('_')[0] for i in samples}
renames['annot_transcript_id'] = 'ID'
df = df.rename(columns=renames)
df = df.set_index('ID')


if snakemake.wildcards['library_prep'] == 'R2C2':
    file_acc = snakemake.params['r2c2'].loc[snakemake.wildcards['sample']].file_acc

    for i in file_acc:
        df[','.join(i)] = df[i].sum(axis=1)
        del df[i[0]]
        del df[i[1]]

df = df / (df.sum(axis=0) / 10**6)

df.to_csv(snakemake.output['expression'], sep='\t')
