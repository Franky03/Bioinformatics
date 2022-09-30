nucleotides= ['A', 'C', 'G', 'T']

def validator(sequence):
    sequence= sequence.upper()
    for n in sequence:
        if n not in nucleotides:
            return False
    return True

def dna_sequence(data):
    if validator(data):
        count_nuc= {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for n in data:
            count_nuc[n] += 1
        return count_nuc
    return False
