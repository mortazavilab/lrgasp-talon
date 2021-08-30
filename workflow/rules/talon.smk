#  Transcript Isoform detection with TALON

rule talon_label_reads:
    input:
        fasta = tc_fasta,
        sam = config['transcriptclean']['sam']
    threads: 16
    resources:
        mem_mb = 64000
    params:
        prefix = config['talon']['sam_read_ann'].replace('_labeled.sam', ''),
        tmp_dir = config['talon']['tmp_dir_read_ann']
    output:
        sam = config['talon']['sam_read_ann']
    shell:
        "talon_label_reads \
        --f {input.sam} \
        --g {input.fasta} \
        --t {threads} \
        --ar 20 \
        --tmpDir={params.tmp_dir} \
        --deleteTmp \
        --o {params.prefix}"


def sam_read_label(wildcards):
    df = read_data_matrix(config['encode']['data_matrix'])

    df = df[
        (df['species'] == wildcards['specie']) &
        (df['sample'] == wildcards['sample']) &
        (df['library_prep'] == wildcards['library_prep']) &
        (df['platform'] == wildcards['platform'])
    ]

    return expand(config['talon']['sam_read_ann'],
                  encode_id=df.index, method=wildcards['method'])


rule talon_config_file:
    input:
        sams = sam_read_label
    output:
        config['talon']['samples']
    threads: 1
    resources:
        mem_mb = 8000
    run:
        name = f"{wildcards.specie}_{wildcards.sample}_{wildcards.library_prep}_{wildcards.platform}"
        df = pd.DataFrame({
            0: name,
            1: [Path(i).stem.replace('_labeled', '') for i in input.sams],
            2: f'{wildcards.library_prep}_{wildcards.platform}',
            3: input.sams
        }).to_csv(output[0], index=False, header=False)


rule talon_initialize_database:
    input:
        gtf = gtf,
    params:
        annot_name = config['talon']['config']['annot_name'],
        genome_name = config['talon']['config']['genome_name'],
        min_transcript_length = config['talon']['config']['min_transcript_length'],
        idprefix = config['talon']['config']['idprefix'],
        dist_5p = config['talon']['config']['dist_5p'],
        dist_3p = config['talon']['config']['dist_3p'],
        db_prefix = config['talon']['db'].replace('.db', '')
    output:
        db = config['talon']['db']
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        "talon_initialize_database \
        --f {input.gtf} \
        --g {params.genome_name} \
        --a {params.annot_name} \
        --l {params.min_transcript_length} \
        --idprefix {params.idprefix} \
        --5p {params.dist_5p} \
        --3p {params.dist_3p} \
        --o {params.db_prefix}"


rule talon_populate_db:
    input:
        config = config['talon']['samples'],
        db = config['talon']['db']
    params:
        tmp_dir = config['talon']['tmp_dir'],
        genome_name = config['talon']['config']['genome_name'],
        prefix = config['talon']['read_annot'].replace(
            '_talon_read_annot.tsv', '')
    output:
        read_annot = config['talon']['read_annot']
    threads: 32
    resources:
        mem_mb = 32000
    shell:
        "talon \
        --f {input.config} \
        --db {input.db} \
        --build {params.genome_name} \
        --tmpDir {params.tmp_dir} \
        -t {threads} \
        --o {params.prefix}"


rule talon_filter_transcripts:
    input:
        db = config['talon']['db'],
        read_annot = config['talon']['read_annot']
    params:
        maxFracA = 0.5,
        minCount = 5,
        minDatasets = 1,
        annot_name = config['talon']['config']['annot_name']
    output:
        white_list = config['talon']['white_list']
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        "talon_filter_transcripts \
        --db {input.db} \
        -a {params.annot_name} \
        --maxFracA={params.maxFracA} \
        --minCount={params.minCount} \
        --minDatasets={params.minDatasets} \
        --o {output.white_list}"


rule talon_gtf:
    input:
        db = config['talon']['db'],
        white_list = config['talon']['white_list']
    params:
        genome_name = config['talon']['config']['genome_name'],
        annot_name = config['talon']['config']['annot_name'],
        gtf_prefix = config['talon']['gtf'].replace('_talon.gtf', '')
    output:
        gtf = config['talon']['gtf']
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        "talon_create_GTF \
        --db {input.db} \
        -b {params.genome_name} \
        -a {params.annot_name} \
        --whitelist={input.white_list} \
        --o {params.gtf_prefix}"
