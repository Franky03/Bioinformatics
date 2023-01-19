import os

file = '3og7.pdb'
file_path = f'{os.getcwd()}/input/{file}'
output_path = f'{os.getcwd()}/hydrogenated/{file}'
reduce = 'C:/Users/kaiky/OneDrive/Documentos/GitHub/Bioinformatics/MolProbity-master/bin/linux/reduce'
cmd = f'{reduce} {file_path} > {output_path}'
os.system(cmd)
