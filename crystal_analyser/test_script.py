# Test very basic functionality of crystal analyser while developing
# Will load a mol2 file of complete molecules generated using the CSD Python API for ABEGEU, having just the molecules in a single unit cell
from Mol2Reader import Mol2Reader
from datetime import datetime
import pandas as pd

def main():
    # Loading in a mol2 file
    start_time = datetime.now()
    input_directory = 'test_script_dir/inputs'
    output_directory = 'test_script_dir/outputs'
    mol2_path = f'{input_directory}/ABEGEU_supercell.mol2'
    reader = Mol2Reader(mol2_path)
    supercell = reader.create_supercell(complete_molecules=True)
    supercell.add_rings_as_atoms()
    # Generating geometric and atomic interactions
    print(f'Number of molecules in the supercell: {len(supercell.molecules)}')
    geometric_interactions = supercell.get_geometric_interactions()
    geometric_interactions.to_csv(f'{output_directory}/geometric_interactions.csv')
    atomic_distances = supercell.get_atomic_distances()
    atomic_distances.to_csv(f'{output_directory}/atomic_distances.csv')
    atomic_interactions = supercell.get_atomic_interactions()
    atomic_interactions.to_csv(f'{output_directory}/atomic_interactions.csv')
    # Getting unique n_mers
    
    # Adding atomic interactions to the supercell

    # Creating graph representations of the supercell

    # 
    print(f'Time taken: {datetime.now() - start_time}')
    
if __name__ == '__main__':
    main()