# Some geometry tools for working with coordinates. Note this is very basic it will be advised to use specfic geometric classes such as Plane
# when i have coded them
import numpy as np

def bond_angle(atom1,atom2,atom3):
    a = atom1.atomic_coordinates
    b = atom2.atomic_coordinates
    c = atom3.atomic_coordinates

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

def vector_angle(v1,v2):
    theta = np.arccos((v1.dot(v2))/(np.sqrt(v1.dot(v1))*np.sqrt(v2.dot(v2))))
    return np.degrees(theta)