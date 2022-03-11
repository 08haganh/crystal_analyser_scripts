# CifReader class
# Only as good as openbabel and pymatgen allow it to be
# Only useful for generating supercells of single crystals based on how Mol2Reader class works
# Better off using other software to generate mol2 files with guarenteed complete molecules
# Pymatgen doesn't necessarily give complete molecules when generating supercell
# Openbabel not always able to identify bonds
# Preserves the labelling of the CIF file by default
import re
import numpy as np
from openbabel import openbabel
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ

class CifReader(CifParser):

    def __init__(self,filename,occupancy_tolerance=1,site_tolerance=0.0001):
        super().__init__(filename,occupancy_tolerance,site_tolerance)
        self.identifier = filename.split('/')[-1].split('.')[0]
        self.cif_dict = self.as_dict()[self.identifier]
        
    def supercell_to_mol2(self,fname,supercell_size,preserve_labelling=True):
        # Can have issues with not writing any bonds to mol2 file
        # however this does not occur often
        name = fname.split('.')[0]
        struc = self.get_structures()[0]
        struc.make_supercell(supercell_size, to_unit_cell=False)
        labels = self.get_new_labels(struc,supercell_size)
        xyzrep = XYZ(struc)
        xyzrep.write_file(f"{name}.xyz")  # write supercell to file
        # Convert supercell to Mol2 format
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("xyz", "mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f"{name}.xyz")   # Open Babel will uncompress automatically
        mol.AddHydrogens()
        obConversion.WriteFile(mol, f'{name}.mol2')
        if preserve_labelling:
            self.change_mol2_atom_labels(f'{name}.mol2',labels)
            
    def supercell_to_xyz(self,fname,supercell_size):
        name = fname.split('.')[0]
        struc = self.get_structures()[0]
        struc.make_supercell(supercell_size, to_unit_cell=False)
        xyzrep = XYZ(struc)
        xyzrep.write_file(f"{name}.xyz")  # write supercell to file
        
    def change_mol2_atom_labels(self,filename,new_labels):
        old_file = open(filename,'r').readlines()
        new_file = open(filename,'w')
        atoms = False
        i=0
        for line in old_file:
            stripped = re.sub("\s+", ",", line.strip())
            split = stripped.split(',')
            arr = np.array(split)
            if arr[0] == '@<TRIPOS>ATOM':
                atoms = True
                new_file.write('@<TRIPOS>ATOM\n')
                continue
            if arr[0] == '@<TRIPOS>BOND':
                atoms = False
                new_file.write('@<TRIPOS>BOND\n')
                continue
            if atoms:
                new_arr = arr
                new_arr[1] = new_labels[i]
                i+=1
            else:
                new_arr = arr
            for elem in new_arr:
                new_file.write(f'{elem} ')
            new_file.write('\n')
        new_file.close()
        
    def get_new_labels(self,struc,supercell_size):
        atom_counter = {}
        new_labels = []
        site_dict = struc.as_dict()['sites']
        symops_len = len(self.cif_dict['_symmetry_equiv_pos_site_id'])
        sc_len = supercell_size[0][0]*supercell_size[1][1]*supercell_size[2][2]
        multiplier = symops_len*sc_len
        for i in range(0,int(len(site_dict)/multiplier)):
            label = site_dict[i*multiplier]['label']
            if label not in atom_counter.keys():
                atom_counter[label] = 1
            new_labels.append([f'{label}{atom_counter[label]}']*multiplier)
            atom_counter[label] += 1
        
        return np.array(new_labels).reshape(-1)