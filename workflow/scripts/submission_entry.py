import json
# import pandas as pd
# from src.utils.io import read_data_matrix


# df = read_data_matrix(snakemake.input['data_matrix'])
# _df = df[
#     (df['species'] == wildcards['specie']) &
#     (df['sample'] == wildcards['sample']) &
#     (df['library_prep'] == wildcards['library_prep']) &
#     (df['platform'] == wildcards['platform'])
# ]

wildcards = snakemake.wildcards

if wildcards['method'] == 'long':
    data_category = 'long_only'
elif wildcards['method'] == 'short':
    data_category = 'long_short'

# remove simulated
submission = {
    "entry_id": f'iso_detect_ref_{wildcards.tool}_{wildcards.method}',
    "challenge_id": snakemake.params['challange_id'],
    "team_name": "3422951",
    "data_category": data_category,
    "samples": [
        "WTC11",
        "H1_mix",
        "ES",
        "mouse_simulation",
        "human_simulation"
    ],
    "library_preps": [
        "cDNA", "dRNA", "R2C2", "CapTrap"
    ],
    "platforms": [
        "ONT", 'PacBio'
    ],
    "experiment_ids": [
        "human_H1mix_CapTrap_PacBio",
        "human_H1mix_cDNA_PacBio",
        "human_H1mix_R2C2_ONT",
        "human_WTC11_CapTrap_PacBio",
        "human_WTC11_cDNA_PacBio",
        "human_WTC11_R2C2_ONT",
        "mouse_ES_CapTrap_PacBio",
        "mouse_ES_cDNA_PacBio",
        "mouse_ES_R2C2_ONT",
        "human_H1mix_CapTrap_ONT",
        "human_H1mix_cDNA_ONT",
        "human_H1mix_dRNA_ONT",
        "human_WTC11_CapTrap_ONT",
        "human_WTC11_cDNA_ONT",
        "human_WTC11_dRNA_ONT",
        "mouse_ES_CapTrap_ONT",
        "mouse_ES_cDNA_ONT",
        "mouse_ES_dRNA_ONT"
    ],
    "contacts": [
        {
            "name": "Muhammed Hasan Celik",
            "email": "mcelik@uci.edu"
        }
    ]
}

with open(snakemake.output['entry'], 'w') as f:
    json.dump(submission, f)
