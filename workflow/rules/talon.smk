# Transcript Isoform detection with TALON

rule talon_label_reads:
    input:
        fasta = fasta_specie,
        sam = config['minimap']['sam'],
        tmp_dir = config['talon']['tmp_dir']
    threads: 16
    params:
        prefix = config['talon']['sam_read_ann'].replace('_labeled.sam', '')
    output:
        sam = config['talon']['sam_read_ann']
    shell:
        "talon_label_reads --f {input.sam} --g {input.fasta} --t {threads} --ar 20 \
        --tmpDir={input.tmp_dir} --deleteTmp --o {params.prefix}"


rule talon_initialize_database:
    input:
        gtf = gtf_specie,
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
        genome_name = config['talon']['config']['genome_name'],
        prefix = config['talon']['read_annot'].replace(
            '_talon_read_annot.tsv', '')
    output:
        read_annot = config['talon']['read_annot']
    threads: 64
    shell:
        "talon \
        --f {input.config} \
        --db {input.db} \
        --build {params.genome_name} \
        --t {threads} \
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
    shell:
        "talon_create_GTF \
        --db {input.db} \
        -b {params.genome_name} \
        -a {params.annot_name} \
        --whitelist={input.white_list} \
        --o {params.gtf_prefix}"
