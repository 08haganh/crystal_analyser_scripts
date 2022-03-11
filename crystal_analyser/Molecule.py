from Bond import Bond
import networkx as nx
import numpy as np
from Plane import Plane
import rdkit.Chem as Chem
from CONFIG import rdkit_bond_types
import pandas as pd
from Ring import Ring
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


class Molecule():

    '''
    The molecule class is the main 
    '''

    def __init__(self,atoms,bonds,canonicalise_atom_order=True):
        self.atoms = atoms
        self.bonds = bonds
        self.has_centre_of_geometry = False
        self.has_plane = False
        self.has_rings = False
        self.plane = self.add_plane()
        self.centre_of_geometry = self.add_centre_of_geometry()
        self.rings = self.add_rings()
        if canonicalise_atom_order:
            self.canonicalise_atom_order()
              
    ############################################### Cool stuff ###########################################
    def add_rings(self):
        self.rings = []
        self.ring_atoms = nx.cycle_basis(self.to_networkx())
        self.ring_bonds = []
        for ring in self.ring_atoms:
            temp = []
            for bond in self.bonds:
                if np.sum(np.isin(bond.atoms,ring)) == 2:
                    temp.append(bond)
                else:
                    continue 
            self.ring_bonds.append(temp)
        for i, (ring_atoms, ring_bonds) in enumerate(zip(self.ring_atoms,self.ring_bonds)):
            for atom in ring_atoms:
                atom.in_ring = True
            for bond in ring_bonds:
                bond.in_ring = True
            ring = Ring(f'ring{i}',ring_atoms,ring_bonds)
            self.rings.append(ring)
        self.has_rings = True
            
    def add_centre_of_geometry(self):
        self.centre_of_geometry = np.average([atom.atomic_coordinates for atom in self.atoms],axis=0)
        self.has_centre_of_geometry = True
    
    def add_plane(self):
        self.plane = Plane(np.array([atom.atomic_coordinates for atom in self.atoms]))
        self.has_plane = True
        
    def get_ring_systems(self):
        ring_systems = Molecule([atom for ring in self.ring_atoms for atom in ring],
                                          [bond for ring in self.ring_bonds for bond in ring],add_rings_as_atoms=False)
        return ring_systems
    
    def get_peripheries(self):
        peripheries = Molecule([atom for atom in self.atoms if (not atom.in_ring)],
                                    [bond for bond in self.bonds if (not bond.in_ring)])
        return peripheries
            
    def test_planarity(self):
        mol_plane = Plane(np.array(atom.atomic_coordinates for atom in self.atoms))
        devs = [mol_plane.point_distance(atom) for atom in self.atoms]
        if np.mean(devs) > 1:
            return False
        else:
            return True

    def get_components(self):
        g = self.to_networkx()
        subgraphs = [g.subgraph(c) for c in nx.connected_components(g)]
        components = []
        for graph in subgraphs:
            bonds = []
            for edge in graph.edges:
                bonds.append(Bond(edge[0],edge[1],g[edge[0]][edge[1]]['type']))
            mol = Molecule(list(graph.nodes),bonds)
            components.append(mol)
        self.components = components
        
        return self.components
        
    def get_unique_components(self):
        g = self.to_networkx()
        subgraphs = [g.subgraph(c) for c in nx.connected_components(g)]
        unique = []
        for i, graph in enumerate(subgraphs):
            if i == 0:
                unique.append(graph)
                continue
            else:
                for un in unique:
                    if nx.isomorphic(un, graph):
                        continue
                    else:
                        unique.append(graph)
        return unique
    
    def canonicalise_atom_order(self):
        atom_labels = np.array([atom.atom_label for atom in self.atoms])
        order = np.argsort(atom_labels)
        self.atoms = np.array(self.atoms)[order].tolist()

    def add_minimum_volume_ellipsoids(self):
        pass
    ############################################### Boring IO stuff ###########################################
    def to_edgelist(self):
        # atom1,atom2,edge_attribute
        pass
    
    def to_bond_dataframe(self):
        # bond,atom1,atom2
        bond_dataframe = []
        for bond in self.bonds:
            bond_dataframe.append({'bond':bond,'atom1':bond.atom1,'atom2':bond.atom2})
        return pd.DataFrame(bond_dataframe)
    
    def to_mol2(self):
        pass
    
    def to_xyz(self,fname):
        split = fname.split('.')
        name = split[0] if len(split) == 1 else split[:-1]
        file = open(name+'.xyz','w')
        n_atoms = len([atom for atom in self.atoms])
        file.write(f'{n_atoms}\n')
        for atom in self.atoms:
            x, y, z = atom.coordinates
            if 'ring' in atom.atomic_symbol:
                file.write(f'Ti {x} {y} {z}\n')
            else:
                file.write(f'{atom.atomic_symbol} {x} {y} {z}\n')
        file.close()

    def to_rdkit(self):
        # adapted from https://github.com/maxhodak/keras-molecules
        mol = Chem.RWMol()
        node_to_idx = {}
        for atom in self.atoms:
            a = Chem.Atom(atom.number)
            idx = mol.AddAtom(a)
            node_to_idx[atom] = idx
        for bond in self.bonds:
            ifirst = node_to_idx[bond.atom1]
            isecond = node_to_idx[bond.atom2]
            bond_type = rdkit_bond_types[bond.bond_type]
            mol.AddBond(ifirst,isecond,bond_type)
        Chem.SanitizeMol(mol)
        return mol.GetMol()

    def to_SMARTS(self):
        mol = self.to_rdkit()
        return Chem.rdmolfiles.MolToSmarts(mol)
    
    def to_networkx(self):
        G = nx.Graph()
        G.add_nodes_from(self.atoms)
        G.add_edges_from([(bond.atom1,bond.atom2,{'type':bond.bond_type}) for bond in self.bonds])
        return G

    def plot(self,ax=None):
        g = self.to_networkx()
        if ax is None:
            fig = plt.figure(figsize=[8,5],dpi=100)
            ax = fig.add_subplot(111,projection="3d")
        else:
            pass
        colour_map = {}
        counter = 0
        for node in g.nodes:
            if node.atomic_symbol not in colour_map.keys():
                colour_map[node.atomic_symbol] = counter
                counter += 1
            x, y, z = node.atomic_coordinates
            ax.scatter(x,y,z,label=node.atom_label,color=cm.tab20([colour_map[node.atomic_symbol]]))
        for edge in g.edges:
            sx, sy, sz = edge[0].atomic_coordinates
            dx, dy, dz = edge[1].atomic_coordinates
            ax.plot((sx,dx),(sy,dy),(sz,dz),label=node.atom_label,color='grey')
        return ax
        

