import argparse
from Bio import SeqIO

# Create argument parser                                                                                                                                                          
parser = argparse.ArgumentParser(prog='Trim sequences to first sequence', usage='%(prog)s -i [aligned_fasta_in]')

# Positional mandatory arguments                                                                                                                                                  
parser.add_argument('-i', '--fasta_in', help = 'Aligned fasta file. First sequence must be ref to trim rest of sequences to', type=str)
args = parser.parse_args()

fasta_in = args.fasta_in
fasta_out = fasta_in.replace('.fasta', '.trimmed.fasta')

fasta_out = fasta_in.replace('.fasta', '.trimmed.fasta')
seqs = list(SeqIO.parse(fasta_in, format = 'fasta'))

ref = seqs[0]
seqs = seqs[1:]

lstrip = len(ref.seq) - len(ref.seq.lstrip('-'))
rstrip = len(ref.seq) - len(ref.seq.rstrip('-'))

ref.seq = ref.seq[lstrip:(len(ref.seq)-rstrip)]

trimmed_seqs_out = [ref]
for s in seqs:
    trimmed_seq = s.seq[lstrip:(len(s.seq)-rstrip)]
    s.seq = trimmed_seq
    trimmed_seqs_out.append(s)
    

SeqIO.write(trimmed_seqs_out, fasta_out, format = 'fasta')

quit('no')
