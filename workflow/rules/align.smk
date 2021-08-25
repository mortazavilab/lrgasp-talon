# Align sort read sequence


def _fq(wildcards):
    df = read_data_matrix(config['encode']['data_matrix'])
    return df[
        (df['species'] == wildcards['specie']) &
        (df['sample'] == wildcards['sample']) &
        (df['platform'] == 'Illumina')
    ]


def fq1(wildcards):
    df = _fq(wildcards)
    return expand(config['encode']['fastq'], encode_id=df.index)


def fq2(wildcards):
    df = _fq(wildcards)
    return expand(config['encode']['fastq'], encode_id=df['paired_acc'])


rule star_index:
    input:
        fasta = fasta,
        gtf = gtf
    output:
        directory(config['star']['index'])
    log:
        "logs/star_index_{specie}.log"
    threads: 32
    resources:
        mem_mb = 128000
    wrapper:
        "0.74.0/bio/star/index"


rule star_align:
    input:
        fq1 = fq1,
        fq2 = fq2,
        index = config['star']['index']
    output:
        bam = config['star']['bam']
    log:
        "logs/star/illumina_{specie}_{sample}_{platform}.log"
    params:
        index = config['star']['index'],
        extra = "--outSAMtype BAM SortedByCoordinate"
    threads: 16
    resources:
        mem_mb = 64000
    wrapper:
        "0.74.0/bio/star/align"


rule gtf2bed:
    input:
        gtf = gtf
    output:
        bed = config['lrgasp']['annotation_bed']
    shell:
        "paftools.js gff2bed {input.gtf} > {output.bed}"


def annotation_bed(wildcards):
    row = _get_row(wildcards['encode_id'])
    return config['lrgasp']['annotation_bed'].format(specie=row['species'])


# Align long read sequencing
rule minimap:
    input:
        ref_fasta = fasta,
        bed = annotation_bed,
        fastq = config['encode']['fastq']
    output:
        sam = config['minimap']['sam']
    log:
        "logs/minimap/{encode_id}.log"
    threads: 16
    resources:
        mem_mb = 64000
    run:
        df = read_data_matrix(config['encode']['data_matrix'])
        row = df.loc[wildcards.encode_id]

        if row.platform == 'PacBio':
            shell("minimap2 --MD -t {threads} -ax splice:hq --junc-bed {input.bed} -uf \
            {input.ref_fasta} {input.fastq} > {output.sam} 2> {log}")
        elif row.platform == 'ONT':
            if row.library_prep == 'dRNA':
                shell("minimap2 --MD -t {threads} -ax splice --junc-bed {input.bed} -k14 -uf \
                {input.ref_fasta} {input.fastq} > {output.sam} 2> {log}")
            else:
                shell("minimap2 --MD -t {threads} -ax splice --junc-bed {input.bed} -uf \
                {input.ref_fasta} {input.fastq} > {output.sam} 2> {log}")


rule sam_to_bam:
    input:
        sam = config['minimap']['sam']
    output:
        bam = config['minimap']['bam']
    threads: 16
    resources:
        mem_mb = 32000
    run:
        shell('samtools view -u {input.sam} |'
              ' samtools sort -@ {threads} -o {output.bam}')
        shell('samtools index {output.bam}')


rule sam_sorted:
    input:
        bam = config['minimap']['bam']
    output:
        sorted_sam = config['minimap']['sam_sorted']
    threads: 1
    resources:
        mem_mb = 8000
    shell:
        "samtools view -h -F 256 -F 4 {input.bam} > {output.sorted_sam}"
