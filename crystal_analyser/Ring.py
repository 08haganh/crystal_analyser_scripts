from RingCentroid import RingCentroid
import numpy as np

class Ring():

    def __init__(self,label,atoms,bonds):
        self.atoms = atoms
        self.bonds = bonds
        self.label = label
        self.aromatic = self.check_aromaticity()
        self.symbol = 'aromatic_ring' if self.aromatic else 'aliphatic_ring'
        self.type = 'aromatic' if self.aromatic else 'aliphatic'
        self.carbon_only = self.check_carbon_only()

    def to_atom(self):
        return RingCentroid(self.atoms,self.label,self.symbol,self.type)

    def check_aromaticity(self):
        lengths = [bond.length() for bond in self.bonds]
        if np.average(lengths) > 1.45: # find paper where you got this number from
            return False
        else:
            return True

    def check_carbon_only(self):
        carbon_only = True
        for atom in self.atoms:
            if atom.atomic_symbol != 'C':
                carbon_only = False
        return carbon_only