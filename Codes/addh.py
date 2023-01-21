import mdtraj as md
from Bio.PDB import PDBParser

pdb = md.load_pdb('./Codes/3og7.pdb')

dssp = md.compute_dssp(pdb)
cont=0
for i, residue in enumerate(pdb.topology.residues):
    if dssp[0][i] == 'NA':
        cont +=1
    print(f"Resíduo {residue} tem estado estrutural secundário {dssp[0][i]}")
print(dssp)
print(cont)