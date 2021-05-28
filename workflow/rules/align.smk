# Align sort read sequence
rule star_index:
    input:
        fasta = fasta_specie,
        gtf = gtf_specie
    output:
        directory(config['star']['index'])
    threads: 16
    log:
        "logs/star_index_{specie}.log"
    wrapper:
        "0.74.0/bio/star/index"


rule star_align:
    input:
        fq1 = config['lrgasp']['simulation']['illumina_fq1'],
        fq2 = config['lrgasp']['simulation']['illumina_fq2'],
        index = config['star']['index']
    output:
        bam = config['star']['bam']
    log:
        "logs/star/illumina_{specie}.log"
    params:
        index = config['star']['index'],
        extra = "--alignSJoverhangMin 8  --alignSJDBoverhangMin 1 --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --twopassMode Basic"
    threads: 16
    wrapper:
        "0.74.0/bio/star/align"


# Align long read sequencing
rule minimap:
    input:
        ref_fasta = fasta_specie,
        fastq = fastq_specie_platform
    output:
        sam = config['minimap']['sam']
    log:
        "logs/minimap/{specie}_{platform}.log"
    threads: 16
    shell:
        "minimap2 --MD -t {threads} -ax splice -uf --secondary=no -C5 \
          {input.ref_fasta} {input.fastq} > {output.sam} 2> {log}"


rule sam_to_bam:
    input:
        sam = config['minimap']['sam']
    output:
        bam = config['minimap']['bam']
    threads: 15  # additional threads
    run:
        shell('samtools view -u {input.sam} |'
              ' samtools sort -@ {threads} -o {output.bam}')
        shell('samtools index {output.bam}')
