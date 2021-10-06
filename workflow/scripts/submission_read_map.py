import pyranges as pr
import pandas as pd


df = pd.read_csv(snakemake.input['read_annot'], sep='\t').rename(
    columns={
        'read_name': 'read_id',
        'annot_transcript_id': 'transcript_id'
    })[['read_id', 'transcript_id']]

df_gtf = pr.read_gtf(snakemake.input['gtf']).df
transcript_ids = df_gtf['transcript_id']
transcript_ids = transcript_ids[~transcript_ids.isna()]
transcript_ids = set(transcript_ids)

df = df[df['transcript_id'].isin(transcript_ids)]
df.to_csv(
    snakemake.output['read_map'], sep='\t',
    index=False, compression='gzip')
