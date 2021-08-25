# Download and prepare row data
rule download_synapse:
    input:
        token = 'configs/synapse_auth_token.txt'
    params:
        project_id = "syn25007493"
    output:
        lrgasp_dir = directory(config['lrgasp']['dir'])
    script:
        "./scripts/synapse_download.py"


rule ungzip_synapse:
    input:
        lrgasp_dir = config['lrgasp']['dir']
    output:
        config['lrgasp']['dir'] + '/../lrgasp.done'
    shell:
        "gzip -rd {input.lrgasp_dir}/* && touch {output}"


rule download_utilities:
    params:
        utilities_dir = config['utilities']['dir'],
        links = [
            "https://raw.githubusercontent.com/LRGASP/lrgasp-submissions/fran/bin/sqanti3_evaluation/utilities/polyA_list.txt",
            "https://raw.githubusercontent.com/LRGASP/lrgasp-submissions/fran/bin/sqanti3_evaluation/utilities/refTSS.mouse.bed",
            "https://raw.githubusercontent.com/LRGASP/lrgasp-submissions/fran/bin/sqanti3_evaluation/utilities/refTSS.human.bed"
        ]
    output:
        config['utilities']['polyA_list'],
        expand(config['utilities']['tss'], specie=config['species'])
    run:
        for i in params['links']:
            shell(f'wget {i} -P {params.utilities_dir}')


rule download_encode:
    input:
        data_matrix = config['encode']['data_matrix']
    output:
        config['encode']['fastq']
    threads: 1
    resources:
        mem_mb = 16000
    script:
        "../scripts/encode_download.py"
