import pdb
import pyranges as pr


gr = pr.read_gtf(snakemake.input['gtf'])
gr = gr[gr.Strand != '.']
gr.to_gtf(snakemake.output['gtf'])
