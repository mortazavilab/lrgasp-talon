rule fasta_header:
    input:
        fasta = fasta,
    output:
        fasta = config['transcriptclean']['fasta']
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/fasta_chrom.py"


def tc_fasta(wildcards):
    wildcards = dict(wildcards)
    if 'specie' in wildcards:
        return config['transcriptclean']['fasta'].format(specie=wildcards['specie'])
    elif 'encode_id' in wildcards:
        row = _get_row(wildcards['encode_id'])
        return config['transcriptclean']['fasta'].format(specie=row['species'])
    else:
        raise ValueError('fasta cannot be determined from wildcards.')


rule transcriptclean_SJ:
    input:
        genome = tc_fasta,
        gtf = gtf
    params:
        tc_path = config['transcriptclean']['tc_path']
    output:
        genome_sj = config['transcriptclean']['gtf_sj']
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        "python {params.tc_path}/accessory_scripts/get_SJs_from_gtf.py \
        --f={input.gtf} \
        --g={input.genome} \
        --o={output.genome_sj}"


def transcriptclean_sj(wildcards):
    row = _get_row(wildcards.encode_id)
    if wildcards.method == 'short':
        return config['star']['SJ'].format(
            specie=row['species'],
            sample=row['sample'],
            platform='Illumina'
        )
    elif wildcards.method == 'long':
        return config['transcriptclean']['gtf_sj'].format(specie=row['species'])
    else:
        raise ValueError('SJ file cannot be determined')


rule transcriptclean_batch:
    input:
        sam = config['minimap']['sam_sorted']
    params:
        batch_size = 100000,
        batch_dir = config['transcriptclean']['batch_dir']
    threads: 1
    resources:
        mem_mb = 4000
    output:
        dynamic(config['transcriptclean']['bam_batch'])
    shell:
        "gatk SplitSamByNumberOfReads \
        -I {input.sam} \
        -O {params.batch_dir} \
        -N_READS {params.batch_size}"


rule transcriptclean_batch_sam:
    input:
        bam = config['transcriptclean']['bam_batch']
    output:
        sam = config['transcriptclean']['sam_batch']
    group: "transcriptclean_batch"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        "samtools view -h {input.bam} > {output.sam}"


rule transcriptclean:
    input:
        sam = config['transcriptclean']['sam_batch'],
        # config['minimap']['sam_sorted'],
        sj = transcriptclean_sj,
        genome = tc_fasta
    params:
        tc_path = config['transcriptclean']['tc_path'],
        tmpdir = config['transcriptclean']['tc_tmpdir'],
        prefix = config['transcriptclean']['sam_clean_batch'].replace(
            '_clean.sam', '')
    output:
        "data/processed/transcriptclean/{encode_id}_{method}/{encode_id}_{method}_batch/{batch_id}/shard_{batch_id}_clean.fa",
        "data/processed/transcriptclean/{encode_id}_{method}/{encode_id}_{method}_batch/{batch_id}/shard_{batch_id}_clean.log",
        "data/processed/transcriptclean/{encode_id}_{method}/{encode_id}_{method}_batch/{batch_id}/shard_{batch_id}_clean.TE.log",
        sam = config['transcriptclean']['sam_clean_batch']
    group: "transcriptclean_batch"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        "python {params.tc_path}/TranscriptClean.py \
        --sam {input.sam} \
        --genome {input.genome} \
        -t {threads} \
        --tmpDir {params.tmpdir} \
        --canonOnly \
        --deleteTmp \
        --spliceJns {input.sj} \
        --outprefix {params.prefix}"


rule transcriptclean_merge:
    input:
        batches = dynamic(config['transcriptclean']['sam_clean_batch'])
    output:
        sam = config['transcriptclean']['sam']
    threads: 4
    resources:
        mem_mb = 32000
    run:
        import os
        tmpdir = f'{resources.tmpdir}/{wildcards.encode_id}_{wildcards.method}'
        if os.path.exists(tmpdir):
            shell('rm -rf {tmpdir}')
        os.mkdir(tmpdir)

        batches_txt = f'{tmpdir}/batches.txt'
        with open(batches_txt, 'w') as f:
            for line in input.batches:
                f.write(line + '\n')
        shell(f"samtools merge - -b {batches_txt} --no-PG | samtools sort -@ {threads} -T {tmpdir} | samtools view -h > {output.sam}")
        shell(f'rm -rf {tmpdir}')
