from Bio.PDB import PDBParser, Selection, NeighborSearch, HSExposureCA
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import calc_angle
import numpy as np
import pandas as pd
import time

#CALCULAR O DSSP ???

class Nodes:
    def __init__(self, name_=None,file_=False):
        self.name= name_
        self.parser= PDBParser(PERMISSIVE=False)
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
        self.ns= NeighborSearch(list(self.structure.get_atoms()))
        

    def Iac(self):

        # Eles não calculam mais essas ligações na nova versão 

        lig_032 = []
        for residue in self.structure.get_residues():
            if str(residue.resname) == "032":
                lig_032.append(residue)
        
        for residue in lig_032:
            for atom in residue:
                for neighbor_pair in self.ns.search(atom.coord, 6.5, level= 'R'):
                    for atom2 in neighbor_pair:
                        if atom2.get_name() == 'CA':
                            distance = np.linalg.norm(atom.coord - atom2.coord)
                            #verificando se o átomo vizinho é de outro resíduo
                            if neighbor_pair != residue:
                                print(residue.resname, neighbor_pair.resname, neighbor_pair.id[1], distance)
    
    def Hidrogen_Bond(self):
        #achar como calcular o angulo entre os átomos
        cutoff = 8.0
        hbond = 3.5
        for residue in self.structure.get_residues():
            for atom in residue:
                neighbors= self.ns.search(atom.coord, cutoff)
                for neighbor in neighbors:
                    distance= np.linalg.norm(atom.coord - neighbor.coord)
                    if atom.fullname[1] in ['N', 'O'] and neighbor.fullname[1] in ['N', 'O']:

                        # terceiro_vetor = ????

                        # terceiro_vetor= carbono_alfa.get_vector()
                        # neighbor_vector= neighbor.get_vector()
                        # a_vector = atom.get_vector()

                        # angle = np.degrees(calc_angle(terceiro_vetor, neighbor_vector, a_vector))

                        if 0.0 < distance <= 3.5:
                            print(atom, residue.id[1], neighbor, neighbor.get_parent().id[1], f"{distance:.3f}")


edges= Edges('file_30', './Codes/3og7.pdb')
edges.Hidrogen_Bond()

