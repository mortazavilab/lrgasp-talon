rule submission_task1_experiment:
    input:
        gtf = config['tools']['gtf'],
        data_matrix = config['encode']['data_matrix']
    output:
        json = config['submission']['task1']['experiment']
    params:
        challange_id = 'iso_detect_ref'
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_structure_task1.py"


rule submission_task2_experiment:
    input:
        gtf = config['tools']['gtf'],
        data_matrix = config['encode']['data_matrix']
    output:
        json = config['submission']['task2']['experiment']
    params:
        challange_id = 'iso_quant'
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_structure_task1.py"


rule benchmark_task1_entry:
    input:
        data_matrix = config['encode']['data_matrix']
    output:
        entry = config['submission']['task1']['entry']
    params:
        challange_id = 'iso_detect_ref'
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_entry.py"

rule benchmark_task2_entry:
    input:
        data_matrix = config['encode']['data_matrix']
    output:
        entry = config['submission']['task2']['entry']
    params:
        challange_id = 'iso_quant'
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_entry.py"


rule benchmark_task1_gtf:
    input:
        gtf = config['tools']['gtf']
    output:
        gtf = config['submission']['task1']['gtf']
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        "cat {input.gtf} | gzip > {output.gtf}"


rule benchmark_task1_read_map:
    input:
        read_annot = config['talon']['read_annot'],
        gtf = config['submission']['task1']['gtf']
    output:
        read_map = config['submission']['task1']['read_map']
    threads: 1
    resources:
        mem_mb = 16000
    script:
        "../scripts/submission_read_map.py"


rule benchmark_task1:
    input:
        submitted_gtf = config['tools']['gtf'],
        gtf = gtf,
        fasta = fasta,
        polyA_list = config['utilities']['polyA_list'],
        tss = config['utilities']['tss'],
        SJ = config['star']['SJ'].format(
            platform='Illumina', specie='{specie}', sample='{sample}'),
        cupcake_path = config['cupcake_path'],
        sqanti3_lrgasp_path = config['sqanti3_lrgasp_path'],
        json = config['submission']['task1']['experiment']
    params:
        name = '{tool}_{specie}_{sample}_{library_prep}_{platform}_{method}',
        conda_env = 'sqanti3_lrgasp'
    output:
        results = directory(config['benchmark_task1']['result_dir'])
    threads: 1
    resources:
        mem_mb = 16000
    run:
        # activate conda env for sqanti3_lrgasp
        # add cupcake and sqanti3_lrgasp to python path
        # run benchmark script
        shell(
            'set +eu '
            ' && PS1=dummy '
            ' && . $(conda info --base)/etc/profile.d/conda.sh '
            ' && conda activate {params.conda_env} '
            ' && echo $CONDA_PREFIX '
            ' && export PYTHONPATH=$PYTHONPATH:{input.cupcake_path}/sequence/ '
            ' && export PYTHONPATH=$PYTHONPATH:{input.cupcake_path} '
            ' && python {input.sqanti3_lrgasp_path} {input.submitted_gtf} {input.gtf} {input.fasta} \
                --gtf --json {input.json} --cage_peak {input.tss} \
                --polyA_motif_list {input.polyA_list} -c {input.SJ} \
                -d {output.results} -o {params.name}'
        )


rule submission_task2_expression:
    input:
        abundance = config['talon']['abundance']
    output:
        expression = config['submission']['task2']['expression']
    params:
        r2c2 = df_r2c2
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_expression.py"


rule benchmark_task2:
    input:
        submitted_gtf = config['tools']['gtf'],
        expression = config['submission']['task2']['expression'],
        evaluation_path = config['benchmark_task2']['path']
    params:
        conda_env = 'lrgasp-challenge-2-evaluation'
    output:
        results = directory(config['benchmark_task2']['result_dir'])
    threads: 8
    resources:
        mem_mb = 32000
    run:
        shell(
            'set +eu '
            ' && PS1=dummy '
            ' && . $(conda info --base)/etc/profile.d/conda.sh '
            ' && conda activate {params.conda_env} '
            ' && echo $CONDA_PREFIX '
            ' && python {input.evaluation_path} \
            -a {input.submitted_gtf} \
            -r {input.expression} \
            -o {output.results} \
            --num_method Single \
            --num_samples Multi'
        )


# def igv_tracks(wildcards):
#     df = read_data_matrix(config['encode']['data_matrix'])
#     df = df[df['sample'] == wildcards['sample']]

#     tracks = [gtf(wildcards)]

#     for i in ['PacBio', 'ONT']:
#         tracks = [*tracks, *expand(
#             config['tools']['gtf'],
#             sample=wildcards['sample'],
#             method=['long', 'short'],
#             library_prep=df['library_prep'].unique(),
#             specie=wildcards['specie'],
#             tool=['talon', 'talon_lapa'],
#             platform=i
#         )]

#     return tracks


# rule igv_report:
#     input:
#         fasta = fasta,
#         tracks = igv_tracks
#     output:
#         "data/processed/igv/{specie}_{sample}.html"
#     log:
#         "logs/igv/{specie}_{sample}_igv-report.log"
#     wrapper:
#         "0.78.0/bio/igv-reports"
