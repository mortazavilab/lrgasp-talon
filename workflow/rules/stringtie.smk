# Transcript Isoform detection with STRINGTIE
rule stringtie_human:
    input:
        bam = config['minimap']['bam'],
        gtf = gtf_specie
    output:
        gtf = config['stringtie']['gtf_raw']
    threads: 4
    shell:
        "stringtie -L -p {threads} -G {input.gtf} -o {output.gtf} {input.bam}"


rule post_process_stringtie:
    input:
        gtf = config['stringtie']['gtf_raw']
    output:
        gtf = config['stringtie']['gtf']
    script:
        "../scripts/post_process_stringtie.py"
