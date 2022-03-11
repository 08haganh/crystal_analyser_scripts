# InteractionDict and Interaction class
import numpy as np
from CONFIG import CONFIG
from GeometryTools import bond_angle, vector_angle

class Interaction():
    def __init__(self,atom1,atom2,mol1_idx=np.nan,mol2_idx=np.nan,atom1_idx=np.nan,atom2_idx=np.nan):
        self.atom1 = atom1
        self.atom2 = atom2
        self.mol1_idx = mol1_idx
        self.mol2_idx = mol2_idx
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.displacement = self.atom2.atomic_coordinates - self.atom1.atomic_coordinates
        self.distance = np.sqrt(np.dot(self.displacement,self.displacement))
        self.vdw_sum = self.atom2.vdw_radius + self.atom1.vdw_radius
        self.vdw_distance = self.distance - self.vdw_sum
        if self.vdw_distance <= 0:
            self.vdw_contact = True
        else:
            self.vdw_contact = False
        self.angle = np.nan
        self.theta1 = np.nan
        self.theta2 = np.nan
        self.vertical_offset = np.nan
        self.horizontal_offset = np.nan
        self.hydrogen_bond_type = np.nan
        self.halogen_bond_type = np.nan
        self.hydrogen_bond = self.check_hydrogen_bond()
        self.halogen_bond = self.check_halogen_bond()
        self.pi_bond = self.check_pi_bond()
        self.ch_pi_bond = self.check_ch_pi_bond()
        self.hydrophobic = self.check_hydrophobic()
        
    def check_hydrogen_bond(self):
        case1 = self.atom1.interaction_dict.hydrogen_bond_donor & self.atom2.interaction_dict.hydrogen_bond_acceptor
        case2 = self.atom2.interaction_dict.hydrogen_bond_donor & self.atom1.interaction_dict.hydrogen_bond_acceptor
        within_distance = ((self.distance < CONFIG['HYDROGEN_BOND']['MAX_DISTANCE']) & 
                            (self.distance > CONFIG['HYDROGEN_BOND']['MIN_DISTANCE']))
        if case1 & within_distance:
            neighbour = self.atom1.neighbours[0]
            angle = bond_angle(neighbour,self.atom1,self.atom2)
            neigh_symbol = neighbour.atomic_symbol 
            if ((angle > CONFIG['HYDROGEN_BOND']['MIN_ANGLE']) & 
                (angle < CONFIG['HYDROGEN_BOND']['MAX_ANGLE'])):
                self.angle = angle
                self.hydrogen_bond_type = neigh_symbol
                return True
            else:
                return False
        elif case2 & within_distance:
            neighbour = self.atom2.neighbours[0]
            angle = bond_angle(neighbour,self.atom2,self.atom1)
            neigh_symbol = neighbour.atomic_symbol 
            if ((angle > CONFIG['HYDROGEN_BOND']['MIN_ANGLE']) & 
                (angle < CONFIG['HYDROGEN_BOND']['MAX_ANGLE'])):
                self.angle = angle
                self.hydrogen_bond_type = neigh_symbol
                return True
            else:
                return False
        else:
            return False
        
    def check_halogen_bond(self):
        # Assign whether halogen bond
        case1 = self.atom1.interaction_dict.halogen_bond_donor & self.atom2.interaction_dict.halogen_bond_acceptor
        case2 = self.atom2.interaction_dict.halogen_bond_donor & self.atom1.interaction_dict.halogen_bond_acceptor
        within_distance = ((self.distance < CONFIG['HALOGEN_BOND']['MAX_DISTANCE']) &
                           (self.distance > CONFIG['HALOGEN_BOND']['MIN_DISTANCE']))
        if (case1 | case2) & within_distance:
            n1 = self.atom1.neighbours[0]
            n2 = self.atom2.neighbours[1]
            theta1 = bond_angle(n1,self.atom1,self.atom2)
            self.theta1 = theta1
            theta2 = bond_angle(n2,self.atom2,self.atom1)
            self.theta2 = theta2
            if ((np.abs(theta2 - theta1) > CONFIG['HALOGEN_BOND']['TYPE1_BOND_DIFFERENCE_MIN']) & 
                (np.abs(theta2 - theta1) < CONFIG['HALOGEN_BOND']['TYPE1_BOND_DIFFERENCE_MAX'])):
                self.halogen_bond_type = 1
            elif ((np.abs(theta2 - theta1) > CONFIG['HALOGEN_BOND']['TYPE1X2_BOND_DIFFERENCE_MIN']) & 
                (np.abs(theta2 - theta1) < CONFIG['HALOGEN_BOND']['TYPE1X2_BOND_DIFFERENCE_MAX'])):
                self.halogen_bond_type = 1.5
            elif ((np.abs(theta2 - theta1) > CONFIG['HALOGEN_BOND']['TYPE2_BOND_DIFFERENCE_MIN']) & 
                (np.abs(theta2 - theta1) < CONFIG['HALOGEN_BOND']['TYPE2_BOND_DIFFERENCE_MAX'])):
                self.halogen_bond_type = 2
            else:
                pass
            return True
        else:
            return False
        
    def check_pi_bond(self):
        # Assign whether pi-pi bond
        case1 = self.atom1.interaction_dict.pi_bond_donor & self.atom2.interaction_dict.pi_bond_acceptor
        case2 = self.atom2.interaction_dict.pi_bond_donor & self.atom1.interaction_dict.pi_bond_acceptor
        within_distance = ((self.distance < CONFIG['PIPI_BOND']['MAX_DISTANCE']) &
                              (self.distance > CONFIG['PIPI_BOND']['MIN_DISTANCE']))
        if (case1 | case2) & within_distance:
            # Calculate bond angle
            # Angle between pi-pi bond and plane of ring1
            pi_plane1 = self.atom1.plane
            pi_plane2 = self.atom2.plane
            pi_bond_angle = pi_plane1.plane_angle(pi_plane2)
            # Calculating offset
            disp = self.atom2.atomic_coordinates - self.atom1.atomic_coordinates
            vec_angle = np.radians(vector_angle(disp, np.array([pi_plane1.a,pi_plane1.b,pi_plane1.c])))
            h_offset = np.abs(self.distance*np.sin(vec_angle))
            v_offset = np.abs(self.distance*np.cos(vec_angle))
            if h_offset < CONFIG['PIPI_BOND']['MAX_OFFSET']:
                if pi_bond_angle > 90:
                    pi_bond_angle = 180 - pi_bond_angle
                within_angle = ((pi_bond_angle > CONFIG['PIPI_BOND']['MIN_ANGLE']) & 
                                (pi_bond_angle < CONFIG['PIPI_BOND']['MAX_ANGLE']))
                if within_angle:
                    self.angle = pi_bond_angle
                    self.horizontal_offset = h_offset
                    self.vertical_offset = v_offset
                    return True
        else:
            return False
    
    def check_ch_pi_bond(self):
        # Assign whether CH-pi bond
        case1 = self.atom1.interaction_dict.ch_pi_bond_donor & self.atom2.interaction_dict.ch_pi_bond_acceptor
        case2 = self.atom2.interaction_dict.ch_pi_bond_donor & self.atom1.interaction_dict.ch_pi_bond_acceptor
        within_distance = ((self.distance < CONFIG['CHPI_BOND']['MAX_DISTANCE']) & 
                           (self.distance > CONFIG['CHPI_BOND']['MIN_DISTANCE']))
        if case1 & within_distance:
            pi_plane = self.atom2.plane
            pi_norm = np.array([pi_plane.a,pi_plane.b,pi_plane.c])
            disp = self.atom2.atomic_coordinates - self.atom1.atomic_coordinates
            pi_bond_angle = np.degrees(np.arccos(disp.dot(pi_norm)/(np.sqrt(disp.dot(disp))*np.sqrt(pi_norm.dot(pi_norm)))))
            if pi_bond_angle > 90:
                pi_bond_angle = 180 - pi_bond_angle
            pi_within_angle = ((pi_bond_angle > CONFIG['CHPI_BOND']['MIN_ANGLE']) & (pi_bond_angle < CONFIG['CHPI_BOND']['MAX_ANGLE']))
            if pi_within_angle:
                self.angle = pi_bond_angle
                return True
        elif case2 & within_distance:
            pi_plane = self.atom1.plane
            pi_norm = np.array([pi_plane.a,pi_plane.b,pi_plane.c])
            disp = self.atom2.atomic_coordinates - self.atom1.atomic_coordinates
            pi_bond_angle = np.degrees(np.arccos(disp.dot(pi_norm)/(np.sqrt(disp.dot(disp))*np.sqrt(pi_norm.dot(pi_norm)))))
            if pi_bond_angle > 90:
                pi_bond_angle = 180 - pi_bond_angle
            pi_within_angle = ((pi_bond_angle > CONFIG['CHPI_BOND']['MIN_ANGLE']) & (pi_bond_angle < CONFIG['CHPI_BOND']['MAX_ANGLE']))
            if pi_within_angle:
                self.angle = pi_bond_angle
                return True
        else:
            return False
        
    def check_hydrophobic(self):
        # Hydrophobic Interactions
        case1 = self.atom1.interaction_dict.hydrophobic & self.atom2.interaction_dict.hydrophobic
        case2 = case1
        within_distance = ((self.distance < CONFIG['CC_HYDROPHOBIC_BOND']['MAX_DISTANCE']) & 
                                       (self.distance > CONFIG['CHPI_BOND']['MIN_DISTANCE']))
        if (case1 | case2) & within_distance:
            return True
        else:
            return False
        
    def to_dict(self):
        info = {
            'mol1_index':self.mol1_idx,
            'mol2_index':self.mol2_idx,
            'atom1_index':self.atom1_idx,
            'atom2_index':self.atom2_idx,
            'atom1_label':self.atom1.atom_label,
            'atom2_label':self.atom2.atom_label,
            'atom1_atomic_symbol':self.atom1.atomic_symbol,
            'atom2_atomic_symbol':self.atom2.atomic_symbol,
            'atom1_type':self.atom1.atom_type,
            'atom2_type':self.atom2.atom_type,
            'a':self.displacement[0],
            'b':self.displacement[1],
            'c':self.displacement[2],
            'distance':self.distance,
            'vdw_sum':self.vdw_sum,
            'vdw_distance':self.vdw_distance,
            'vdw_contact':self.vdw_contact,
            'hydrogen_bond':self.hydrogen_bond,
            'halogen_bond':self.halogen_bond,
            'pi_bond':self.pi_bond,
            'ch_pi_bond':self.ch_pi_bond,
            'hydrophobic':self.hydrophobic,
            'angle':self.angle,
            'theta1':self.theta1,
            'theta2':self.theta2,
            'horizontal_offset':self.horizontal_offset,
            'vertical_offset':self.vertical_offset,
            'hydrogen_bond_type':self.hydrogen_bond_type,
            'halogen_bond_type':self.halogen_bond_type}
        
        return info    