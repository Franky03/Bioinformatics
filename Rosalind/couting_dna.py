def dna_sequence(data):
    count_nuc= {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for n in data:
        if n in count_nuc.keys():
            count_nuc[n] += 1
    return count_nuc

with open('./Rosalind/rosalind_dna.txt', 'r') as sequence:
    dna= sequence.read()

print(' '.join([str(val) for key, val in dna_sequence(dna).items()]))