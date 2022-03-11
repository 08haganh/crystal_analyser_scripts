# Atom class
from CONFIG import vdw_radii, atomic_mass, atomic_number
from InteractionDict import InteractionDict

class Atom():
    
    '''
    Atom class created from Mol2Reader
    The role of the Atom class is to be able to map the atomic interactions that are present in a supercell
    To do this, knowledge of the van der Waals radius, atomic symbol, atom type, and atom neighbours are required
    All of this information is recorded upon supercell creation using the Mol2Reader
    '''

    def __init__(self,atom_label,atomic_coordinates,atomic_symbol,atom_type):
        self.atom_label = atom_label # str :: 'C1'
        self.atomic_coordinates = atomic_coordinates # np.array :: np.array([0,0,0],dtype=np.float32)
        self.atomic_symbol = atomic_symbol # str :: 'C'
        self.atom_type = atom_type # str :: 'C.ar
        self.interaction_dict = False # InteractionDict class
        self.in_ring = False # bool :: True / False
        self.neighbours = [] # list of atom objects :: needed to add interaction_dict
        try:
            self.atomic_mass = atomic_mass[self.symbol]
        except:
            self.atomic_mass = 0
        try:
            self.atomic_number = atomic_number[self.symbol]
        except:
            self.atomic_number = 0
        try:
            self.vdw_radius = vdw_radii[self.symbol]
        except:
            self.vdw_radius = 0

    def add_interaction_dict(self):
        self.interaction_dict = InteractionDict(self)