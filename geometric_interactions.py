import os 
import argparse
import numpy as np
from tqdm import tqdm

import sys
sys.path.append('crystal_analyser')
from Mol2Reader import Mol2Reader

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--mol2_directory',type=str)
    parser.add_argument('--n_splits',type=int)
    parser.add_argument('--split_index',type=int)
    parser.add_argument('--output_directory',type=str)
    parser.add_argument('--output_prefix',type=str)
    parser.add_argument('--n_processes',default=4,type=int)
    args = parser.parse_args()
    mol2_directory = args.mol2_directory
    n_splits = args.n_splits
    split_index = args.split_index
    output_directory = args.output_directory
    output_prefix = args.output_prefix
    n_processes = args.n_processes
    try:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
    except:
        pass
    # get mol2 paths
    mol2_paths = os.listdir(mol2_directory)
    mol2_paths.sort()
    # Get split
    splits = np.array_split(np.array(mol2_paths,dtype=str),n_splits)
    split = splits[split_index-1]
    print(split_index-1)
    print(len(splits))
    print(len(split))
    # Loop through mol2 files and calculate geometric interactions
    for path in tqdm(split):
        try:
            identifier = path.split('.')[0]
            reader = Mol2Reader(os.path.join(mol2_directory,path))
            supercell = reader.create_supercell(complete_molecules=True)
            interactions = supercell.get_geometric_interactions(n_processes)
            interactions.to_csv(f'{output_directory}/{output_prefix}{identifier}_geometric_interactions.csv')
        except Exception as e:
            print(f'{identifier} has failed')

if __name__ == '__main__':
    main()