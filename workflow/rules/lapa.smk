
rule chrom_sizes:
    input:
        fasta = fasta
    threads: 1
    resources:
        mem_mb = 8000
    output:
        chrom_sizes = config['lrgasp']['chrom_sizes']
    shell:
        "faidx {input.fasta} -i chromsizes > {output.chrom_sizes}"


def bams(wildcards):
    df = read_data_matrix(config['encode']['data_matrix'])
    df = df[
        (df['species'] == wildcards['specie']) &
        (df['sample'] == wildcards['sample']) &
        (df['library_prep'] == wildcards['library_prep']) &
        (df['platform'] == wildcards['platform'])
    ]
    return expand(config['minimap']['bam'], encode_id=df.index)


rule lapa:
    input:
        bams = bams,
        fasta = fasta,
        gtf = gtf,
        chrom_sizes = chrom_sizes
    threads: 1
    resources:
        mem_mb = 32000
    output:
        lapa_dir = directory(config['lapa']['result_dir'])
    run:
        bam = ','.join(input.bams)
        shell(f"lapa \
        --alignment {bam} \
        --fasta {input.fasta} \
        --annotation {input.gtf} \
        --chrom_sizes {input.chrom_sizes} \
        --method tail \
        --output_dir {output}")


rule lapa_correct_gtf:
    input:
        read_annot = config['talon']['read_annot'],
        lapa_dir = config['lapa']['result_dir'],
        gtf = config['talon']['gtf'],
        fasta = fasta
    threads: 1
    resources:
        mem_mb = 32000
    output:
        gtf = config['lapa']['gtf']
    script:
        "../scripts/correct_gtf.py"
