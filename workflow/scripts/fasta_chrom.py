from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


seqs = list()

for record in tqdm(SeqIO.parse(snakemake.input['fasta'], 'fasta')):
    seqs.append(SeqRecord(record.seq, id=record.id, description=''))

SeqIO.write(seqs, snakemake.output['fasta'], 'fasta')
