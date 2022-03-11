# Example script for basic crystal analyser usage
import os
from CifReader import CifReader
from Atom import *
from Bond import *
from RingCentroid import *
from Molecule import *
from Ring import *
from Supercell import *

def main():
    # Generating a supercell from a cif file
    cif_reader = CifReader('example_files/ABEGEU.cif')
    supercell_size = [[5,0,0],[0,5,0],[0,0,5]]
    cif_reader.supercell_to_mol2('example_files/crystal_analyser_ABEGEU_supercell_555.mol2',supercell_size,preserve_labelling=True)
    # Reading a single crystal mol2 file generated using CrystalAnalyser
    simple_supercell = Mol2Reader('example_files/crystal_analyser_ABEGEU_supercell_555.mol2',complete_molecules=False)
    # Reading a co-crystal mol2 supercell created in the CSD
    co_supercell = Mol2Reader('example_files/csd_ABEGEU_supercell_555.mol2',complete_molecules=True)
    # Adding interactions to single crystal

    # Adding interactions to co-crystal
    # Writing single crystal with interactions to mol2 file
    # Writing co-crystal with interactions to mol2 file
    # Drawing 

if __name__ == '__main__':
    main()
