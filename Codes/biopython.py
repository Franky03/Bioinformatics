from Bio.PDB import PDBParser, Selection, NeighborSearch
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Residue import Residue
import numpy as np
import pandas as pd
import time

#CALCULAR O DSSP ???

class Nodes:
    def __init__(self, name_=None,file_=False):
        self.name= name_
        self.parser= PDBParser()
        self.structure= self.parser.get_structure(name_, file_)
        self.model= self.structure[0]
        #self.dssp = DSSP(self.model, './Codes/3og7.pdb', file_type='PDB', dssp='dssp')
        self.all_dssps= []
        self.nodes_id, self.chains, self.positions, self.residues = [], [], [], []
        #Degrees não está com uma precisão boa, é preciso saber como pode ser feito igual o do RINGs
        self.degrees= []
        self.cut_dist= 8.0 # definir um limite de distância de corte (escolhido com base na literatura)
        #B-Factor, coords and filenames
        self.bfactors, self.coords, self.pdb_filename =[], [], []
        self.models= []
    
    def dssp_test(self):
        self.search_nodes()
        for residue in self.dssp:
            print(residue)
    
    def search_nodes(self):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if str(residue.resname) != 'HOH': # ignore solvent
                        # Node_ID, Chain, Position and Residue
                        self.nodes_id.append(f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}")
                        self.chains.append(str(chain.id))
                        self.positions.append(residue.id[1])
                        if str(residue.resname) == '032':
                            #Se o resíduo for o 032 então ele não terá um Bfactor-CA nem coordenadas
                            self.bfactors.append('NaN')
                            self.coords.append(np.array(['NaN', 'NaN', 'NaN']))
                        self.residues.append(str(residue.resname))
                        
                        #DSSP not working

                        self.all_dssps.append('NaN')
                        #Degree

                        degree= 0

                        if 'CA' in residue:
                        #R ->  list of residues list of modules
                            for residue_2 in Selection.unfold_entities(model, 'R'):
                                if ( residue_2.get_id()[1] != residue.get_id()[1] ):
                                    #CA -> Carbono Alfa (mais preciso que outros átomos)
                                    #Calcular a distancia euclidiana das coordenadas desse carbono alfa
                                    if 'CA' in residue_2:
                                        distance= np.linalg.norm(residue["CA"].coord - residue_2["CA"].coord)
                                        if ( distance < self.cut_dist ):
                                            degree += 1 
                            
                        self.degrees.append(degree)

                        #Bfactor_CA
                        b_factor=0
                        count=0
                        for atom in residue:
                            if ( atom.get_name() == 'CA' ):
                                b_factor += atom.get_bfactor()
                                count += 1

                                coords= atom.get_coord()
                                self.coords.append(coords)
                                
                        if (count != 0):
                            bf_average= b_factor/count
                            self.bfactors.append(f"{bf_average:.3f}")
                        
                        #pdb filenames

                        self.pdb_filename.append(f"input_file.pdb#{str(residue.id[1])}.{str(chain.id)}")
                        self.models.append(model.id + 1)

    def print_output(self):
        start= time.time()
        self.search_nodes()
        for n in range(len(self.nodes_id)):
            try:
                print(f"{self.nodes_id[n]}\t{self.chains[n]}\t\t{self.positions[n]}\t\t{self.residues[n]}\t{self.all_dssps[n]}\t" +
                        f"{self.degrees[n]}\t{self.bfactors[n]:.3f}\t{self.coords[n][0]:.3f}\t{self.coords[n][1]:.3f}\t{self.coords[n][2]:.3f}\t" +
                        f"{self.pdb_filename[n]}\t{self.models[n]}"
                )
            except Exception as e:
                print(f"{self.nodes_id[n]}\t{self.chains[n]}\t\t{self.positions[n]}\t\t{self.residues[n]}\t{self.all_dssps[n]}\t" +
                        f"{self.degrees[n]}\t{self.bfactors[n]}\t{self.coords[n][0]}\t{self.coords[n][1]}\t{self.coords[n][2]}\t" +
                        f"{self.pdb_filename[n]}\t{self.models[n]}"
                )
        print(f"---{(time.time() - start)} seconds ---")

    def to_file(self):
        start= time.time()
        self.search_nodes()

        colunas= ["NodeId", "Chain","Position",	"Residue",	"Dssp",	"Degree", "Bfactor_CA", "x", "y", "z", "pdbFileName", "Model"]
        x = [f"{x_c[0]:.3f}" for x_c in self.coords]
        y = [f"{y_c[1]:.3f}" for y_c in self.coords]
        z = [f"{z_c[2]:.3f}" for z_c in self.coords]

        data= pd.DataFrame(list(zip(self.nodes_id, self.chains, self.positions, self.residues, 
                                    self.all_dssps, self.degrees, self.bfactors, x, y, z, self.pdb_filename, self.models
        )), columns= colunas)
        
        data.to_csv(f'./{self.name}.csv', sep='\t', index=False)

        print(f"---{(time.time() - start)} seconds ---")

class Edges(Nodes):
    def __init__(self, name, file_pdb):
        Nodes.__init__(self ,name_=name, file_= file_pdb)
        self.edges = []
        self.res= [res for res in self.structure.get_residues()]
        self.params = {
            'Hydrogen_Bond': [3.1, 0],
            'Salt_Bridge': [4.0, 0],
            'Dissulfide_Bond': [2.2, 0],
            'Van_der_Waals': [3.2, 0],
            'Pi_Stacking': [7.2, 0],
            'Sulfur_Aryl': [7.2, 0],
            'Cation_Aryl': [6, 0],
            'Anion_Aryl': [3.5, 0],
            "tshaped": 0,
            "inter": 0,
            "paralel": 0,
            'lighbacep': ['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1'],
            'lighbdono': ['HE1', 'H', 'HE', 'HH11', 'HH12', 'HH21',
                          'HH22', 'HD1', 'HE2', 'HG', 'HG1', 'HG21', 'HG22',
                          'HG23', 'HH', 'HD21', 'HD22', 'HE21', 'HE22'],
            'ligsb1': ['OD2', 'OE2', 'OH'],
            'ligsb2': ['NZ', 'NH2', 'NH1'],
            'aasbneg': ['ASP', 'GLU'],
            'aasbpos': ['LYS', 'ARG'],
            'ligvdw': ['CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE'],
            'aavdw': ['VAL', 'TRE', 'MET', 'LEU', 'ILE'],
            'aa_arom': ['TYR', 'PHE', 'TRP'],
            'ligctn': ['MG', 'CU', 'K', 'FE2', 'FE', 'NI', 'NA', 'MO1', 'MO3', 'MO4',
                       'MO5', 'MO6', 'MO7', 'MO8', 'MO9', 'NZ', 'NH2', 'NH1'],
            'ligctn2': ['CG', 'CE2', 'CG'],
            'aactn_beg': 3.4,
            'aactn_end': 4.0,
            'ligan1': ['CL', 'BR', 'I', 'OD2', 'OE2', 'OH'],
            'ligspi1': ['SG'],
            'ligspi2': ['CG', 'CE2', 'CG'],
            'aaan_beg': 2.0,
            'aaan_end': 3.0,
            'aaspi': 5.3
            }
    
    def search_edges(self):
        for i, res1 in enumerate(self.res):
            if 'CA' in res1:
                for j, res2 in enumerate(self.res):
                    if i<j:
                        if 'CA' in res2:
                            distance= np.linalg.norm(res1["CA"].coord - res2["CA"].coord)
                            print(res1.resname, res2.resname)
    
    def bonds(self):
        ns= NeighborSearch(list(self.structure.get_atoms()))
        bonds= ns.search_all(3)
        for bond in bonds:
            atom1, atom2= bond
            print(atom1, atom2)
            if atom1.element == 'H' or atom2 == 'H':
                print('HBOND')
            elif (atom1.element == "C" and atom2.element == "C") or (atom1.element == "N" and atom2.element == "N"):
                print('VDW')

                    

node= Edges('file_30', './Codes/file_30.pdb')
node.bonds()

