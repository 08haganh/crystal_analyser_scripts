from CONFIG import CONFIG
import numpy as np

class InteractionDict():
    def __init__(self,atom):
        self.atom = atom
        self.check_hydrogen_bond_donor()
        self.check_hydrogen_bond_acceptor()
        self.check_halogen_bond_donor()
        self.check_halogen_bond_acceptor()
        self.check_pi_bond_donor()
        self.check_pi_bond_acceptor()
        self.check_ch_pi_bond_donor()
        self.check_ch_pi_bond_acceptor()
        self.check_hydrophobic()
        
    def check_hydrogen_bond_donor(self):
        if self.atom.atomic_symbol == 'H':
            neighbours = [atom.atomic_symbol for atom in self.atom.neighbours]
            assert len(neighbours) > 0
            if  np.sum(np.isin(np.array(neighbours),np.array(CONFIG['HYDROGEN_BOND']['DONORS']))) > 0:
                self.hydrogen_bond_donor = True 
            else:
                self.hydrogen_bond_donor = False
        else:
            self.hydrogen_bond_donor = False
        
    def check_hydrogen_bond_acceptor(self):
        if self.atom.atomic_symbol in CONFIG['HYDROGEN_BOND']['ACCEPTORS']:
            self.hydrogen_bond_acceptor = True 
        else:
            self.hydrogen_bond_acceptor = False
            
    def check_halogen_bond_donor(self):
        if self.atom.atomic_symbol in CONFIG['HALOGEN_BOND']['DONORS']:
            self.halogen_bond_donor = True
        else:
            self.halogen_bond_donor = False
            
    def check_halogen_bond_acceptor(self):
        if self.atom.atomic_symbol in CONFIG['HALOGEN_BOND']['ACCEPTORS']:
            self.halogen_bond_acceptor = True
        else:
            self.halogen_bond_acceptor = False
        
    def check_pi_bond_donor(self):
        if self.atom.atomic_symbol in CONFIG['PIPI_BOND']['DONORS']:
            self.pi_bond_donor = True
        else:
            self.pi_bond_donor = False 
            
    def check_pi_bond_acceptor(self):
        if self.atom.atomic_symbol in CONFIG['PIPI_BOND']['ACCEPTORS']:
            self.pi_bond_acceptor = True
        else:
            self.pi_bond_acceptor = False 
            
    def check_ch_pi_bond_donor(self):
        if self.atom.atomic_symbol in CONFIG['CHPI_BOND']['DONORS']:
            neighbours = neighbours = [atom.atomic_symbol for atom in self.atom.neighbours]
            assert len(neighbours) > 0
            if  np.sum(np.isin(np.array(neighbours),np.array(['C']))) > 0:
                self.ch_pi_bond_donor = True
            else:
                self.ch_pi_bond_donor = False
        else:
            self.ch_pi_bond_donor = False
    
    def check_ch_pi_bond_acceptor(self):
        if self.atom.atomic_symbol in CONFIG['CHPI_BOND']['ACCEPTORS']:
            self.ch_pi_bond_acceptor = True
        else:
            self.ch_pi_bond_acceptor = False
            
    def check_hydrophobic(self):
        if self.atom.atomic_symbol == 'C':
            neighbours = neighbours = [atom.atomic_symbol for atom in self.atom.neighbours]
            assert len(neighbours) > 0
            if  np.sum(np.isin(np.array(neighbours),np.array(['C','H']),invert=True)) == 0:
                self.hydrophobic = True
            else:
                self.hydrophobic = False
        else:
            self.hydrophobic = False