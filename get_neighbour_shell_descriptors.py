# Script for generating nearest neighbour shells for the n nearest shells from a directory of geometric interactions
# Output will be a single dataframe indexed by the filename with n entries for the n nearest neighbour shells
# The vector describing each shell will be the number of neighbours in the shell followed by the centroid distance, interplanar angle, vertical offset, and horizontal offset
import os 
import argparse
import pandas as pd
import numpy as np

import concurrent.futures

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_directory',type=str)
    parser.add_argument('--n_shells',type=int,default=4)
    parser.add_argument('--output_directory',type=str)
    parser.add_argument('--output_prefix',type=str)
    parser.add_argument('--n_processes',type=int,default=6)
    args = parser.parse_args()
    input_directory = args.input_directory
    n_shells = args.n_shells
    output_directory = args.output_directory
    output_prefix = args.output_prefix
    n_processes = args.n_processes
    # Create output directory
    try:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
    except:
        pass
    # Load file paths
    filenames = os.listdir(input_directory)
    filenames.sort()
    # Create batches for multiprocessing
    filenames = np.array(filenames,dtype=object).reshape(-1,1)
    input_directory_arr = np.array([input_directory]*len(filenames),dtype=object).reshape(-1,1)
    n_shells_arr = np.array([n_shells]*len(filenames),dtype=object).reshape(-1,1)
    batches = np.concatenate([filenames,input_directory_arr,n_shells_arr],axis=1)
    # Loop through files and get nearest neighbour shell dataframes
    print(f'Getting neighbour shell descriptors from {len(filenames)} geometric interactions')
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes) as ex:
        results = ex.map(get_neighbour_shells,batches)
    pd.concat(results).to_csv(f'{output_directory}/{output_prefix}_nshells{n_shells}_shell_descriptors.csv')

def get_neighbour_shells(batch):
    fname = batch[0]
    input_directory = batch[1]
    n_shells = batch[2]
    columns = ['centroid_distance','interplanar_angle','vertical_offset','horizontal_offset']
    geometric_interactions = pd.read_csv(os.path.join(input_directory,fname),index_col=[0,1])[columns]
    geometric_interactions = geometric_interactions.loc[geometric_interactions.centroid_distance < 20]
    geometric_interactions = sanitise_interactions(geometric_interactions)
    geometric_interactions = geometric_interactions.sort_values('centroid_distance',ascending=True)
    unique_interactions = geometric_interactions.drop_duplicates(keep='first')
    geometric_interactions = geometric_interactions.reset_index().set_index(columns)
    unique_interactions = unique_interactions.reset_index().set_index(columns)
    neighbour_shells = []
    for i, unique_interaction in enumerate(unique_interactions.index):
        if i == n_shells:
            break
        temp_inter = np.array(unique_interaction).reshape(len(columns),1)
        mol1_counts = geometric_interactions['mol1_index'].sort_index().loc[unique_interaction].value_counts()
        mol2_counts = geometric_interactions['mol2_index'].sort_index().loc[unique_interaction].value_counts()
        counts = pd.DataFrame(mol1_counts).join(pd.DataFrame(mol2_counts)).sum(axis=1)
        max_neighbours = np.array(max(counts)).reshape(-1,1)
        desc = np.concatenate([max_neighbours,temp_inter]).reshape(len(columns)+1)
        neighbour_shells.append(desc)

    return pd.DataFrame(neighbour_shells,index=[fname]*n_shells)

def sanitise_interactions(interactions):
        # Removes small numerical differences between interactions by replacing with first incident
        # required for matching interactions for joins, unique interactions, etc
        length = len(interactions.columns)
        seen = np.array([np.zeros(shape=(length))])
        new_values = []
        mask = []
        for idx in interactions.index:
            values = interactions.loc[idx].values
            if list(values) in seen.tolist():
                mask.append(False)
            else:
                if np.sum((np.sum(np.isclose(values,seen,atol=0.05),axis=1) == length),axis=0)>0:
                    mask.append(True)
                    new_values.append(seen[np.sum(np.isclose(values,seen,atol=0.05),axis=1) == length][0])   
                else:
                    mask.append(False)
                    seen = np.append(seen,[values], axis=0)
        interactions[mask] = new_values
        return interactions

if __name__ == '__main__':
    main()