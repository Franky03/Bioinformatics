from BioTools import *
import random

DNA_STR= ''.join([random.choice(nucleotides) for n in range(30)])

print(dna_sequence(DNA_STR))