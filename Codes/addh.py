import pymol

pymol.cmd.load('3og7.pdb', 'myprotein')
pymol.cmd.h_add()
pymol.cmd.save('./hydrogenated/30g7.pdb')
