rule submission_task1_experiment:
    input:
        gtf = config['tools']['gtf'],
        data_matrix = config['encode']['data_matrix']
    output:
        json = config['submission']['task1']['experiment']
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/submission_structure_task1.py"


rule benchmark_task1_entry:
    input:
        data_matrix = config['encode']['data_matrix']
    output:
        config['submission']['task1']['entry']
    threads: 1
    resources:
        mem_mb = 4000
    run:
        "../scripts/submission_task1_entry.py"


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
