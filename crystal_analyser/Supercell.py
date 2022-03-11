import numpy as np
import pandas as pd
import concurrent.futures
from Interaction import Interaction
from GeometryTools import vector_angle

N_PROCESSES = 4

class Supercell(object):

    def __init__(self,molecules):
        self.molecules = molecules

    ################################################################################### IO ###############################################################################################
    def to_centroid_graph(self):
        pass

    def to_mol2(self,add_interactions=False):
        pass

    ############################################################################### STRUCTURE ###############################################################################################
    def add_rings_as_atoms(self):
        for mol in self.molecules:
            temp_atoms = []
            for ring in mol.rings:
                ring_atom = ring.to_atom()
                ring_atom.add_interaction_dict()
                temp_atoms.append(ring_atom)
            mol.atoms = mol.atoms + temp_atoms

    def remove_rings_as_atoms(self):
        for mol in self.molecules:
            mol.atoms = [atom for atom in mol.atoms if 'ring' not in atom.symbol]

    ############################################################################## INTERACTIONS ###############################################################################################
    def get_geometric_interactions(self,N_PROCESSES=N_PROCESSES):
        # Create ordered list of mol1, mol2, mol1_index, mol2_index for multiprocessing self.geometric_interaction()
        mol1s = []
        mol2s = []
        mol1_indices = []
        mol2_indices = []
        for i, mol1 in enumerate(self.molecules,0):
            for j, mol2 in enumerate(self.molecules[i+1:],i+1):
                mol1s.append(mol1)
                mol2s.append(mol2)
                mol1_indices.append(i)
                mol2_indices.append(j)
        mol1s = np.array(mol1s).reshape(-1,1)
        mol2s = np.array(mol2s).reshape(-1,1)
        mol1_indices = np.array(mol1_indices).reshape(-1,1)
        mol2_indices = np.array(mol2_indices).reshape(-1,1)
        batches = np.concatenate([mol1s,mol2s,mol1_indices,mol2_indices],axis=1)
        batches = np.array_split(batches,N_PROCESSES)
        print(f'Generating geometric interactions for {batches[0].shape[0]} by {len(batches)} for {len(batches)*batches[0].shape[0]} pairs of molecules')
        with concurrent.futures.ProcessPoolExecutor(max_workers=N_PROCESSES) as executor:
            interactions = executor.map(self.geometric_interaction,batches)
        return np.round(pd.DataFrame([item for sublist in interactions for item in sublist]).set_index(['mol1_index','mol2_index']),6)

    def geometric_interaction(self,batch):
        interactions = []
        for entry in batch:
            mol1 = entry[0]
            mol2 = entry[1]
            mol1_index = entry[2]
            mol2_index = entry[3]
            # Calculate entry interaction
            disp = mol2.centre_of_geometry - mol1.centre_of_geometry
            dist = np.sqrt(np.dot(disp,disp))
            interplanar_angle = mol2.plane.plane_angle(mol1.plane)
            vec_angle = np.radians(vector_angle(disp, mol1.plane.unit_normal()))
            v_offset = np.abs(dist*np.cos(vec_angle))
            h_offset = np.abs(dist*np.sin(vec_angle))
            interaction = {'mol1_index':mol1_index,
                        'mol2_index':mol2_index,
                        'centroid_x':disp[0],
                        'centroid_y':disp[1],
                        'centroid_z':disp[2],
                        'centroid_distance':dist,
                        'interplanar_angle':interplanar_angle,
                        'vertical_offset':v_offset,
                        'horizontal_offset':h_offset}
            interactions.append(interaction)

        return interactions


    def get_atomic_interactions(self,cutoff=5.5):
        atomic_distances = self.get_atomic_distances().reset_index()
        atomic_distances = atomic_distances.loc[atomic_distances.distance <= cutoff]
        batches = np.array_split(atomic_distances.values,N_PROCESSES)
        print(f'Generating atomic interactions for {batches[0].shape[0]} by {len(batches)} for {len(batches)*batches[0].shape[0]} pairs of atoms')
        with concurrent.futures.ProcessPoolExecutor(max_workers=N_PROCESSES) as executor:
            interactions = executor.map(self.atomic_interaction,batches)
        return np.round(pd.DataFrame([item for sublist in interactions for item in sublist]).set_index(['mol1_index','mol2_index']),6)

    def atomic_interaction(self,batch):
        interactions = []
        for entry in batch:
            mol1_index = int(entry[0])
            mol2_index = int(entry[1])
            atom1_index = int(entry[2])
            atom2_index = int(entry[3])
            atom1 = self.molecules[mol1_index].atoms[atom1_index]
            atom2 = self.molecules[mol2_index].atoms[atom2_index]
            interactions.append(Interaction(atom1,atom2,mol1_index,mol2_index,atom1_index,atom2_index).to_dict())

        return interactions

    def get_atomic_distances(self):
        mol1s = []
        mol2s = []
        mol1_indices = []
        mol2_indices = []
        for i, mol1 in enumerate(self.molecules,0):
            for j, mol2 in enumerate(self.molecules[i+1:],i+1):
                mol1s.append(mol1)
                mol2s.append(mol2)
                mol1_indices.append(i)
                mol2_indices.append(j)
        mol1s = np.array(mol1s).reshape(-1,1)
        mol2s = np.array(mol2s).reshape(-1,1)
        mol1_indices = np.array(mol1_indices).reshape(-1,1)
        mol2_indices = np.array(mol2_indices).reshape(-1,1)
        batches = np.concatenate([mol1s,mol2s,mol1_indices,mol2_indices],axis=1)
        batches = np.array_split(batches,N_PROCESSES)
        print(f'Calculating atomic distances for {batches[0].shape[0]} by {len(batches)} for {len(batches)*batches[0].shape[0]} pairs of molecules')
        columns = ['mol1_index','mol2_index','atom1_index','atom2_index','distance']
        with concurrent.futures.ProcessPoolExecutor(max_workers=N_PROCESSES) as executor:
            interactions = executor.map(self.atomic_distance,batches)
        return np.round(pd.DataFrame((np.vstack([item for sublist in interactions for item in sublist])),
                                    columns=columns).set_index(['mol1_index','mol2_index']),6)

    def atomic_distance(self,batch):
        # Calculate atomic distances between pair of molecules on batch
        # pads numpy arrays if molecules have different number of atoms
        # faster than nested for loops
        interactions = []
        for entry in batch:
            mol1 = entry[0]
            mol2 = entry[1]
            mol1_index = entry[2]
            mol2_index = entry[3]
            atom1_indices = []
            atom2_indices = []
            atom1_coordinates = []
            atom2_coordinates = []
            for i, atom in enumerate(mol1.atoms):
                atom1_indices.append(i)
                x,y,z = atom.atomic_coordinates
                atom1_coordinates.append((x,y,z))
            for i, atom in enumerate(mol2.atoms):
                atom2_indices.append(i)
                x,y,z = atom.atomic_coordinates
                atom2_coordinates.append((x,y,z))
            if len(atom1_coordinates) >= len(atom2_coordinates):
                small_atoms = len(atom2_coordinates)
                small_mol_index = np.array([mol2_index]*small_atoms)
                small_atom_indices = np.array(atom2_indices)
                small_atom_coordinates = np.array(atom2_coordinates)
                big_mol_index = np.array([mol1_index]*small_atoms)
                big_atoms = len(atom1_coordinates)
                big_atom_indices = np.array(atom1_indices)
                big_atom_coordinates = np.array(atom1_coordinates)
            else:
                small_atoms = len(atom1_coordinates)
                small_mol_index = np.array([mol1_index]*small_atoms)
                small_atom_indices = np.array(atom1_indices)
                small_atom_coordinates = np.array(atom1_coordinates)
                big_mol_index = np.array([mol2_index]*small_atoms)
                big_atoms = len(atom2_coordinates)
                big_atom_indices = np.array(atom2_indices)
                big_atom_coordinates = np.array(atom2_coordinates)
            for x in range(big_atoms):
                padding_coords = np.zeros(shape=(big_atoms-small_atoms,3))
                small_atom_coordinates = np.concatenate([small_atom_coordinates,padding_coords],axis=0)
                disp = big_atom_coordinates - small_atom_coordinates
                dist = np.sqrt(disp[:,0]*disp[:,0] + disp[:,1]*disp[:,1] + disp[:,2]*disp[:,2]).reshape(-1,1)
                dist = dist[:small_atoms]
                big_atom_indices_reduced = big_atom_indices[:small_atoms]
                # Ensure smaller index always in column 1
                if big_mol_index[0] <= small_mol_index[0]:
                    interactions.append(
                        np.concatenate([big_mol_index.reshape(-1,1),
                                        small_mol_index.reshape(-1,1),
                                        big_atom_indices_reduced.reshape(-1,1),
                                        small_atom_indices.reshape(-1,1),
                                        dist],axis=1))
                else:
                    interactions.append(
                        np.concatenate([small_mol_index.reshape(-1,1),
                                        big_mol_index.reshape(-1,1),
                                        small_atom_indices.reshape(-1,1),
                                        big_atom_indices_reduced.reshape(-1,1),
                                        dist],axis=1))
                big_atom_indices = np.roll(big_atom_indices,-1,axis=0)
                big_atom_coordinates = np.roll(big_atom_coordinates,-1,axis=0)
        return np.array(interactions)