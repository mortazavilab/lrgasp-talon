# export TMPDIR='/pub/mcelik/tmp'

# printf "Using temp directory: "
# python -c "import os; print(os.environ['TMPDIR'])"

python -m snakemake \
       -j 500 \
       --keep-going \
       --latency-wait 30 \
       --cluster "sbatch --mem {resources.mem_mb}M -c {threads} -A seyedam_lab --cluster-constraint=fastscratch​3 --time=9:00:00"       
