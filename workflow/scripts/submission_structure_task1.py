import json
from src.utils.io import read_data_matrix

wildcards = snakemake.wildcards

sample = wildcards.sample

if sample == 'H1mix':
    sample = 'H1_mix'

library_prep = [wildcards.library_prep]

if wildcards.method == 'short':
    library_prep.append('cDNA')

platform = [wildcards.platform]

if wildcards.method == 'short':
    platform.append('Illumina')


df = read_data_matrix(snakemake.input['data_matrix'])
_df = df[
    (df['species'] == wildcards['specie']) &
    (df['sample'] == wildcards['sample']) &
    (df['library_prep'] == wildcards['library_prep']) &
    (df['platform'] == wildcards['platform'])
]
libraries = _df.index.tolist()

if wildcards.method == 'short':
    _df = df[
        (df['species'] == wildcards['specie']) &
        (df['sample'] == wildcards['sample']) &
        (df['platform'] == 'Illumina')
    ]
    libraries = [*libraries, *_df.index.tolist(), *_df['paired_acc'].tolist()]
    platform.append('Illumina')

if wildcards.method == 'long':
    data_category = 'long_only'
elif wildcards.method == 'short':
    data_category = 'long_short'

submission = {
    "experiment_id": f'{wildcards.specie}_{wildcards.sample}_{wildcards.library_prep}_{wildcards.platform}',
    "challenge_id": snakemake.params['challange_id'],
    "description": f'{wildcards.tool}_{wildcards.method}_{wildcards.specie}_{wildcards.sample}_{wildcards.library_prep}_{wildcards.platform}',
    "notes": "",
    "species": wildcards.specie,
    "data_category": data_category,
    "samples": [sample],
    "library_preps": list(set(library_prep)),
    "platforms": platform,
    "libraries": libraries,
    "software": [
        {
            "name": "TALON",
            "description": "A technology-agnostic long-read analysis pipeline for transcriptome discovery and quantification",
            "version": "4.0.0",
            "url": "https://github.com/mortazavilab/lrgasp-talon",
            "config": "All parameters are defined workflow in the github link"
        },
        {
            "name": "LAPA",
            "description": "Alternative polyadenylation detection from diverse data sources such as 3'-seq, long-read and short-reads.",
            "version": "0.0.2",
            "url": "https://github.com/mortazavilab/lapa",
            "config": "All parameters are defined workflow in the github link"
        }
    ]
}


with open(snakemake.output['json'], 'w') as f:
    json.dump(submission, f)
