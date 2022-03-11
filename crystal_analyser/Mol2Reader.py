# Class for reading mol2 files
# returns a list of Molecule objects with geometrical features e.g. centres of geometry and molecular planes
# chemical structure information such as ring systems
import re
from Atom import Atom
from Bond import Bond
from Molecule import Molecule
from Supercell import Supercell
import networkx as nx 
import numpy as np

class Mol2Reader():
    '''
    Mol2Reader Class - Reads a mol2 file from which it creates a supercell object
                     - works by constructing a networkx graph from atom and bond information within mol2 file
                     - then for each disconnected subgraph, creates rdkit molecules to check its chemical validity
                     - if chemically valid, then the subgraph is kept, else discarded
    '''
    def __init__(self,path):
        self.path = path
        self.file = open(self.path,'r')

    def create_supercell(self,complete_molecules=False):
        atoms, bonds = self.create_atoms_and_bonds() # Reads mol2 file to generate self.atoms and self.bonds
        molecules = self.create_molecules(atoms,bonds,complete_molecules=complete_molecules) # Creates molecules and rings using networkx
        return Supercell(molecules)

    def create_atoms_and_bonds(self):
        atoms = []
        bonds = []
        tripos_atom = False
        tripos_bond = False
        for line in self.file.readlines():
            arr = self.line_to_array(line)
            if '@<TRIPOS>' in arr[0]:
                if arr[0] == '@<TRIPOS>ATOM':
                    tripos_atom = True
                    continue
                if arr[0] == '@<TRIPOS>BOND':
                    tripos_atom = False
                    tripos_bond = True
                    continue
                else:
                    tripos_atom = False
                    tripos_bond = False
                    continue
            # Create list of Atom objects
            if tripos_atom:
                atom_number = int(arr[0])
                atom_label = str(arr[1])
                x = float(arr[2])
                y = float(arr[3])
                z = float(arr[4])
                atom_type = str(arr[5])
                atom_coordinates = np.array([x,y,z])
                atom_symbol = re.sub("\d+", "",atom_label)
                atoms.append(Atom(atom_label,atom_coordinates,atom_symbol,atom_type))
            # Create list of Bond objects; assign atom neighbours using bonding information
            if tripos_bond:
                bond_atom_number_1 = int(arr[1])
                bond_atom_number_2 = int(arr[2])
                bond_type = str(arr[3])
                bond_atom1 = atoms[bond_atom_number_1-1]
                bond_atom2 = atoms[bond_atom_number_2-1]
                atoms[bond_atom_number_1-1].neighbours.append(atoms[bond_atom_number_2-1])
                atoms[bond_atom_number_2-1].neighbours.append(atoms[bond_atom_number_1-1])
                bonds.append(Bond(bond_atom1,bond_atom2,bond_type))

        return (atoms,bonds)
        
    def create_molecules(self,atoms,bonds,complete_molecules):
        molecules = []
        supergraph = nx.Graph()
        supergraph.add_nodes_from(atoms)
        supergraph.add_edges_from([(bond.atom1,bond.atom2,{'type':bond.bond_type}) for bond in bonds])
        subgraphs = [supergraph.subgraph(c) for c in nx.connected_components(supergraph)]
        # Using n_atoms potentially buggy
        # Will have to have a think as to how to load co-crystals
        n_atoms = max([len(subgraph.nodes) for subgraph in subgraphs])
        if not complete_molecules:
            subgraphs = [subgraph for subgraph in subgraphs if len(subgraph.nodes) == n_atoms]
        else:
            subgraphs = subgraphs
        for graph in subgraphs:
            bonds = []
            for edge in graph.edges:
                bonds.append(Bond(edge[0],edge[1],supergraph[edge[0]][edge[1]]['type']))
            mol = Molecule(list(graph.nodes),bonds)
            molecules.append(mol)
        for mol in molecules:
            mol.add_centre_of_geometry()
            mol.add_plane()
            mol.add_rings()
            for atom in mol.atoms:
                atom.add_interaction_dict()

        return molecules
         
    def line_to_array(self,line):
        stripped = re.sub("\s+", ",", line.strip())
        split = stripped.split(',')
        arr = np.array(split)
        return arr