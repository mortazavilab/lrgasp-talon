import json


df = read_data_matrix(snakemake.input['data_matrix'])
_df = df[
    (df['species'] == wildcards['specie']) &
    (df['sample'] == wildcards['sample']) &
    (df['library_prep'] == wildcards['library_prep']) &
    (df['platform'] == wildcards['platform'])
]

if wildcards.method == 'long':
    data_category = 'long_only'
elif wildcards.method == 'short':
    data_category = 'long_short'


submission = {
    "entry_id": "iso_detect_ref_{wildcards.tool}_{wildcards.method}",
    "challenge_id": "iso_detect_ref",
    "team_id": "syn1",
    "team_name": "MortazaviLab",
    "data_category": data_category,
    "samples": [
        "WTC11",
        "H1_mix",
        "ES",
        "mouse_simulation"
    ],
    "library_preps": [
        "dRNA"
    ],
    "platforms": [
        "ONT"
    ],
    "experiment_ids": [
        "ES_drna_ont",
        "H1_mix_drna_ont",
        "mouse_sim_drna_ont",
        "WTC11_drna_ont"
    ],
    "contacts": [
        {
            "name": "Muhammed Hasan Celik",
            "email": "mcelik@uci.edu"
        }
    ]
}
