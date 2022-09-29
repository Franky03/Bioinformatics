def dna_sequence(data):
    if len(data)== 0:
        return 0

    A= data.count('A')
    C= data.count('C')
    G= data.count('G')
    T = data.count('T')

    return (A, C, G , T)

print(dna_sequence('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'))
        