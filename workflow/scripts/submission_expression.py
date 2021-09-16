import pandas as pd


df = pd.read_csv(snakemake.input['abundance'], sep='\t')
samples = df.columns[11:].tolist()
df = df[['annot_transcript_id'] + samples]

renames = {i: i.split('_')[0] for i in samples}
renames['annot_transcript_id'] = 'ID'
df.rename(columns=renames).to_csv(
    snakemake.output['expression'], sep='\t', index=False)
