
python -m snakemake \
       -j 500 \
       --keep-going \
       --latency-wait 60 \
       --cluster "sbatch --mem {resources.mem_mb}M -c {threads} -A seyedam_lab --cluster-constraint=fastscratch​3"
