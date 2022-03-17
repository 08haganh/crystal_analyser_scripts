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
    #filenames = np.array(['BARWUO_geometric_interactions.csv','BARWUO_geometric_interactions.csv'],dtype=object).reshape(-1,1)
    input_directory_arr = np.array([input_directory]*len(filenames),dtype=object).reshape(-1,1)
    n_shells_arr = np.array([n_shells]*len(filenames),dtype=object).reshape(-1,1)
    batches = np.concatenate([filenames,input_directory_arr,n_shells_arr],axis=1)
    # Uncomment code below when not testing
    #get_neighbour_shells(batches[0]).to_csv(f'{output_directory}/{output_prefix}nshells{n_shells}_shell_descriptors_test.csv')
    # Loop through files and get nearest neighbour shell dataframes
    print(f'Getting neighbour shell descriptors from {len(filenames)} geometric interactions')
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes) as ex:
        results = ex.map(get_neighbour_shells,batches)
    pd.concat(results).to_csv(f'{output_directory}/{output_prefix}nshells{n_shells}_shell_descriptors.csv')

def get_neighbour_shells(batch):
    fname = batch[0]
    input_directory = batch[1]
    n_shells = batch[2]
    return_columns = ['centroid_distance','interplanar_angle','vertical_offset','horizontal_offset','displacement_plane_dp']
    sort_columns = ['centroid_distance','interplanar_angle']#,'vertical_offset','horizontal_offset','displacement_plane_dp']
    geometric_interactions = pd.read_csv(os.path.join(input_directory,fname),index_col=[0,1])
    geometric_interactions = geometric_interactions.loc[geometric_interactions.centroid_distance < 20]
    geometric_interactions = np.round(geometric_interactions,2)
    reduced_geometric_interaction = geometric_interactions[sort_columns].copy()
    geometric_interactions = sanitise_interactions(geometric_interactions)
    geometric_interactions = geometric_interactions.sort_values('centroid_distance',ascending=True)
    unique_interactions = geometric_interactions.drop_duplicates(keep='first')
    reduced_geometric_interaction = geometric_interactions[sort_columns].copy()
    reduced_unique_interactions = reduced_geometric_interaction.drop_duplicates(keep='first')
    reduced_geometric_interaction = reduced_geometric_interaction.reset_index().set_index(sort_columns)
    reduced_unique_interactions = reduced_unique_interactions.reset_index().set_index(sort_columns)
    unique_interactions = unique_interactions.reset_index().set_index(sort_columns)
    neighbour_shells = []
    for i, reduced_unique_interaction in enumerate(reduced_unique_interactions.index):
        if i == n_shells:
            break
        unique_locator = np.array(reduced_unique_interaction).reshape(len(sort_columns))
        #print(unique_locator)
        #print(unique_interactions.sort_index())
        #print(reduced_geometric_interaction['mol1_index'].sort_index())
        temp_inter = unique_interactions.sort_index().loc[reduced_unique_interaction]
        #print(temp_inter.displacement_plane_dp)
        dpdp = np.average(np.unique(temp_inter.displacement_plane_dp))
        if type(temp_inter) == pd.DataFrame:
            temp_inter = temp_inter.iloc[0]
        #print(temp_inter)
        cd, ipa = temp_inter.name
        temp_inter = pd.DataFrame(temp_inter.values,index=temp_inter.index).T
        temp_inter[sort_columns] = [cd, ipa]
        temp_inter = temp_inter[return_columns].values.reshape(len(return_columns))
        mol1_counts = reduced_geometric_interaction['mol1_index'].sort_index().loc[reduced_unique_interaction].value_counts()
        mol2_counts = reduced_geometric_interaction['mol2_index'].sort_index().loc[reduced_unique_interaction].value_counts()
        counts = pd.DataFrame(mol1_counts).join(pd.DataFrame(mol2_counts)).sum(axis=1)
        max_neighbours = np.array(max(counts)).reshape(1)
        #yh1 # print(max_neighbours,'\n',temp_inter)
        desc = np.concatenate([max_neighbours,temp_inter]).reshape(len(return_columns)+1)
        desc[-1] = dpdp
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