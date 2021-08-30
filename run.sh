python -m snakemake \
       -j 500 \
       --keep-going \
       --cluster "sbatch --mem {resources.mem_mb}M -c {threads} -A seyedam_lab"
       # --default-resources "tmpdir=/pub/mcelik/tmpdir"

