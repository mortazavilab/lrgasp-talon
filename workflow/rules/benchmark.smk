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
        sqanti3_lrgasp_path = config['sqanti3_lrgasp_path']
    params:
        name = '{tool}_{specie}_{sample}_{library_prep}_{platform}_{method}',
        conda_env = 'sqanti3_lrgasp'
    output:
        results = directory(config['benchmark']['result_dir'])
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
                --gtf --name {params.name} \
                --platform {wildcards.platform} --cage_peak {input.tss} \
                --polyA_motif_list {input.polyA_list} -c {input.SJ} \
                -d {output.results} -o {params.name}'
        )

# rule benchmark_task2
# "https://github.com/Tidesun/LRGASP_visualization"


# TODO: benchmark task 2
# TODO: benchmark task 3
