import numpy as np

class Bond(object):

    '''
    Bond class 
    Required for assigning atomic neighbours
    '''

    def __init__(self,atom1,atom2,bond_type):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_type = bond_type
        self.atoms = [atom1,atom2]

    def length(self):
        disp = self.atom2.atomic_coordinates - self.atom1.atomic_coordinates
        return np.sqrt(disp.dot(disp))
