from Atom import Atom
import numpy as np
from Plane import Plane

class RingCentroid(Atom):

    def __init__(self,atoms,atom_label,atom_symbol,atom_type):
        coordinates = np.average(np.array([atom.atomic_coordinates for atom in atoms]),axis=0)
        self.plane = Plane(np.array([atom.atomic_coordinates for atom in atoms]))
        super().__init__(atom_label,coordinates,atom_symbol,atom_type)
