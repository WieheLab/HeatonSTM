from Bio import SeqIO
import pandas as pd
import string
import argparse
import time
import re

## script to get AA frequency by position when aligned to a reference. Originally designed to use a reference sequence to align to and another sequence to use as the numbering scheme (but these may be the same sequence)

letters = list(string.ascii_uppercase)

codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'stop', 'TAG':'stop',
    'TGC':'C', 'TGT':'C', 'TGA':'stop', 'TGG':'W'}

# Create argument parser                                                                                                                                                          
parser = argparse.ArgumentParser(prog='AA frequency matrix by position', usage='%(prog)s -i [aligned_fasta_in] -N [name_out] -R [fasta with reference and numbering ref aligned (numbering ref first)]')

# Positional mandatory arguments                                                                                                                                                  
parser.add_argument('-i', '--fasta_in', help = 'Aligned fasta file. First sequence must be ref to trim rest of sequences to', type=str)
parser.add_argument('-N',"--name_out", help="Prefix to name output csv files", type=str, default = '')
parser.add_argument('-R', '--aligned_ref_fasta', help='Fasta with numbering reference and reference aligned (numbering reference listed first), if both are the same then the fasta needs to have the sequences twice', type = str)
args = parser.parse_args()

fasta_in = args.fasta_in
name_out = args.name_out
aligned_ref_fasta = args.aligned_ref_fasta


### Get numbering from fasta with aligned ref and numbering ref
def get_hxb2_numbering(aligned_ref_fasta):
    refs = list(SeqIO.parse(aligned_ref_fasta, format = 'fasta'))

    hxb2 = refs[0]
    ref = refs[1]

    i = 1
    letters_num = 0
    pos_num = []
    hxb2_out = []
    ref_out = []
    for h,r in zip(hxb2.seq, ref.seq):
        i_in = i
        if h != '-':
            if r!='-':
                pos_num.append(i)
                hxb2_out.append(h)
                ref_out.append(r)      
            i += 1
        else:
            num = str(i-1) + letters[letters_num]
            pos_num.append(num)
            hxb2_out.append(h)
            ref_out.append(r)
        if i != i_in:
            letters_num = 0
        else:
            letters_num +=1

    hxb2_numbering_out = pd.DataFrame([range(1,len(ref_out)+1), pos_num, ref_out], index = ['Reference Position', 'Numbered Reference Position', 'Reference AA'])
    print('Got HXB2 numbering...')
    return(hxb2_numbering_out)

hxb2_numbering_out = get_hxb2_numbering(aligned_ref_fasta)
hxb2_numbering_out

# trim to reference and remove sequences that don't overlap with reference
def trim_and_remove_nonoverlapping_seqs(ref, seqs):
    lstrip = len(ref.seq) - len(ref.seq.lstrip('-'))
    rstrip = len(ref.seq) - len(ref.seq.rstrip('-'))

    ref.seq = ref.seq[lstrip:(len(ref.seq)-rstrip)]

    trimmed_seqs_out = []
    for s in seqs:
        trimmed_seq = s.seq[lstrip:(len(s.seq)-rstrip)]
        if re.search('[AGCT]', str(trimmed_seq))!= None:
            s.seq = trimmed_seq
            trimmed_seqs_out.append(s)
    return(trimmed_seqs_out)
    

print('Reading in fasta...')
seqs = list(SeqIO.parse(fasta_in, format = 'fasta'))
print('Finished reading in fasta...')
print(len(seqs), 'sequences found...')

ref = seqs[0]
ref_ungap = str(ref.seq).replace('-','')
seqs = seqs[1:]

seqs = trim_and_remove_nonoverlapping_seqs(ref, seqs)
len(seqs)




def get_aa_codon_column_names(ref_ungap):
    column_names_aa = []
    column_names_codon = []

    i = 0
    k=1
    while i < len(ref_ungap):
        ref_codon = ref_ungap[int(i):int(i+3)]
        ref_aa = codon_table[ref_codon]
        column_names_aa.append(str(k) + ':' + ref_aa)
        column_names_codon.append(str(k) + ':' +ref_codon)
        i += 3
        k += 1
    return(column_names_aa, column_names_codon)

column_names_aa, column_names_codon = get_aa_codon_column_names(ref_ungap)

def get_ungapped_seq_aligned_to_ref(s):
    ungap_seq = []
    
    for ref_char,char in zip(str(ref.seq),str(s.seq)):
        if ref_char != '-':
            ungap_seq.append(char)
            
    seq_codon_list = []
    seq_aa_list = []
    
    i = 0
    while i < len(ref_ungap):
        seq_codon = ungap_seq[int(i):int(i+3)]
        seq_codon = ''.join(seq_codon)
        if '-' in seq_codon:
            seq_aa = '-'
        elif len(set(seq_codon) - set("AGCT")) == 0:
            seq_aa = codon_table[seq_codon]
        seq_codon_list.append(seq_codon)
        seq_aa_list.append(seq_aa) 
        i += 3
    
    return(seq_codon_list, seq_aa_list)

def get_new_codon_list(codon_list_in):
    # get left missing count
    flag = 0
    left_count = 0
    for r in codon_list_in:
        if flag == 0:
            if r == '---':
                left_count +=1
            if r != '---':
                flag = 1 
    #get right missing count
    
    out_codon_reverse = codon_list_in[::-1]
    flag = 0
    right_count = 0
    for r in out_codon_reverse:
        if flag == 0:
            if r == '---':
                right_count +=1
            if r != '---':
                flag = 1 
    
    for r in range(0, left_count):
        codon_list_in[r] = 'MISSING'

    for r in range(1, right_count+1):
        codon_list_in[-r] = 'MISSING'
    return(codon_list_in)

def get_new_aa_list(aa_list_in):
    # get left missing count
    flag = 0
    left_count = 0
    for r in aa_list_in:
        if flag == 0:
            if r == '-':
                left_count +=1
            if r != '-':
                flag = 1 
    #get right missing count
    
    out_aa_reverse = aa_list_in[::-1]
    flag = 0
    right_count = 0
    for r in out_aa_reverse:
        if flag == 0:
            if r == '-':
                right_count +=1
            if r != '-':
                flag = 1 
    
    for r in range(0, left_count):
        aa_list_in[r] = 'MISSING'

    for r in range(1, right_count+1):
        aa_list_in[-r] = 'MISSING'
    return(aa_list_in)




codon_list_out = []
aa_list_out = []
tic = time.perf_counter()
for s in seqs:
    codon_out, aa_out = get_ungapped_seq_aligned_to_ref(s)
    new_codon_out = get_new_codon_list(codon_out)
    new_aa_out = get_new_aa_list(aa_out)
    codon_list_out.append(new_codon_out)
    aa_list_out.append(new_aa_out)
    
toc = time.perf_counter()
print(f"Got codon and AA lists in {toc - tic:0.4f} seconds")

mat_codon = pd.DataFrame(codon_list_out, columns = column_names_codon)
mat_aa = pd.DataFrame(aa_list_out, columns = column_names_aa)




print('Converted codon and AA list to dataframes...')

row_names = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 
           'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'stop', '-', 'MISSING']
mat_aa_freqs = pd.DataFrame(0,columns = column_names_aa, index = row_names)

row_names_codon = codon_table.keys()
mat_codon_freqs = pd.DataFrame(0, columns = column_names_codon, index = row_names_codon)

for col in mat_aa.columns:
    aa_list = list(mat_aa[col])
    freqs_aa = {n: float(aa_list.count(n)) for n in set(aa_list)}
    for x in freqs_aa.keys():
        mat_aa_freqs.loc[x,col] = freqs_aa[x]
    
mat_aa_freqs.drop(axis = 0, index = 'MISSING', inplace = True)
mat_aa_freqs = mat_aa_freqs.fillna(0)
print('Finished AA freqs matrix...')



sums = []
for col in mat_aa_freqs.columns:
    sums.append(mat_aa_freqs[col].sum())
    mat_aa_freqs[col] = mat_aa_freqs[col] / mat_aa_freqs[col].sum()

mat_aa_freqs =  pd.concat([mat_aa_freqs, pd.DataFrame([sums], index = ['Total'], columns= mat_aa_freqs.columns )])

for col in mat_codon.columns:
    codon_list = list(mat_codon[col])
    freqs_codon = {n: float(codon_list.count(n)) for n in set(codon_list)}
    for x in freqs_codon.keys():
        if '-' in x:
            mat_codon_freqs.loc["'"+''.join(x),col] = freqs_codon[''.join(x)]
        else:
            mat_codon_freqs.loc[x,col] = freqs_codon[x]

mat_codon_freqs.drop(axis = 0, index = 'MISSING', inplace = True)
mat_codon_freqs = mat_codon_freqs.fillna(0)
print('Finished codon freqs matrix...')

sums = []
for col in mat_codon_freqs.columns:
    sums.append(mat_codon_freqs[col].sum())
    mat_codon_freqs[col] = mat_codon_freqs[col] / mat_codon_freqs[col].sum()
    
mat_codon_freqs = pd.concat([mat_codon_freqs, pd.DataFrame([sums], index = ['Total'], columns= mat_codon_freqs.columns )])

mat_codon_freqs

print('Writing to csv...')

hxb2_numbering_out.columns = mat_aa_freqs.columns
mat_aa_freqs_out = pd.concat([hxb2_numbering_out, mat_aa_freqs])

mat_aa_freqs_out.to_csv(name_out +'.matrix_out_aa.csv', index = True, header = False)
mat_codon_freqs.to_csv(name_out + '.matrix_out_codons.csv')

pos_list = []
ref_aa_list = []
aa_list = []
freq_list = []
mat_aa_freqs_out.drop(axis = 0, index = 'Total', inplace = True)
col_list = mat_aa_freqs_out.columns
row_list = [x for x in mat_aa_freqs_out.index[3:]]
hxb2_nums = []
for c in col_list:
    pos = c.split(':')[0]
    ref_aa = c.split(':')[1]
    if ref_aa not in ['stop', '-']:
        for r in row_list:
            if r not in ['stop', '-'] and r!= ref_aa:
                aa_list.append(r)
                freq = mat_aa_freqs.loc[r, c]
                freq_list.append(freq)
                ref_aa_list.append(ref_aa)
                pos_list.append(pos)
                hxb2 = mat_aa_freqs_out.loc['Numbered Reference Position', c]
                hxb2_nums.append(hxb2)

            

df_new = pd.DataFrame({'Position': pos_list, 'Numbered Reference Position': hxb2_nums, 'Reference AA': ref_aa_list, 'AA': aa_list, 'Frequency':freq_list})
convert_dict = {'Position': str,
                'Numbered Reference Position': str,
                'Reference AA': str,
                'AA': str,
                'Frequency': float
               }
  
df_new = df_new.astype(convert_dict)
df_new = df_new.sort_values(by = ['Frequency'], axis = 0, ascending = False)

df_new.to_csv(name_out + '.mutation_frequencies.csv', index = False)
            


