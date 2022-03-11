# Script for generating nearest neighbour shells for the n nearest shells from a directory of geometric interactions
# Output will be a single dataframe indexed by the filename with n entries for the n nearest neighbour shells
# The vector describing each shell will be the number of neighbours in the shell followed by the centroid distance, interplanar angle, vertical offset, and horizontal offset
import os 
import argparse
import pandas as pd
import sympy as sym
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
    batches = np.concatenate([filenames,input_directory_arr],axis=1)
    # Loop through files and get nearest neighbour shell dataframes
    print(batches[0])
    get_neighbour_shells(batches[0]).to_csv(f'{output_directory}/{output_prefix}_fundamental_neighbour_shells_test.csv')
    print(f'Getting neighbour shell descriptors from {len(filenames)} geometric interactions')
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes) as ex:
        results = ex.map(get_neighbour_shells,batches)
    pd.concat(results).to_csv(f'{output_directory}/{output_prefix}_fundamental_neighbour_shells.csv')

def get_neighbour_shells(batch):
    fname = batch[0]
    input_directory = batch[1]
    columns = ['centroid_distance','interplanar_angle','vertical_offset','horizontal_offset']
    data = pd.read_csv(os.path.join(input_directory,fname),index_col=[0,1])
    data = np.round(data,3)
    data = data.sort_values('centroid_distance')
    data = data.loc[data.centroid_distance < 20]
    data = sanitise_interactions(data)
    unique = data.drop_duplicates(keep='first')
    centroid_displacements = unique[['centroid_x','centroid_y','centroid_z']]
    centroid_displacements.shape
    distances = data[columns]
    basis = get_independent_basis(centroid_displacements)
    distances = np.round(np.sqrt(basis[:,0]*basis[:,0] + basis[:,1]*basis[:,1] + basis[:,2]*basis[:,2]),3)
    unique_interactions = []
    for d in distances:
        if np.any(np.isclose(unique.centroid_distance,d,atol=0.05)):
            unique_interactions.append(pd.DataFrame(unique.loc[np.isclose(unique.centroid_distance,d,atol=0.05)].iloc[0]).T)
        else:
            pass
    unique_interactions = pd.concat(unique_interactions)
    unique_interactions.drop_duplicates()
    unique_interactions = unique_interactions.drop_duplicates()
    neighbour_shells = []
    for cd in unique_interactions.centroid_distance:
        temp_inter = unique.loc[unique.centroid_distance == cd][columns]
        mol1_counts = data.loc[data.centroid_distance == cd].index.get_level_values(0).value_counts()
        mol2_counts = data.loc[data.centroid_distance == cd].index.get_level_values(1).value_counts()
        counts = pd.DataFrame(mol1_counts).join(pd.DataFrame(mol2_counts)).sum(axis=1)
        desc = np.concatenate([np.array(max(counts)).reshape(-1,1),temp_inter.values[0].reshape(-1,1)],axis=0).reshape(-1)
        neighbour_shells.append(desc)
    
    return pd.DataFrame(neighbour_shells,index=[fname]*len(neighbour_shells))

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

def check_basis_combinations(basis,vector,dim=3):
    vector = vector.reshape(dim,1)
    lock_i = False
    lock_j = False
    lock_k = False
    for ic in range(basis.shape[1]-2):
        i = i if lock_i else ic
        #i = ic
        for jc in range(ic+1,basis.shape[1]-1):
            j = j if lock_j else jc
            #j = jc
            for kc in range(jc+1,basis.shape[1]):
                k = k if lock_k else kc
                #k = kc
                #print(i,j,k)
                temp_basis = np.round(np.concatenate([basis[:,i].reshape(dim,1),basis[:,j].reshape(dim,1),basis[:,k].reshape(dim,1)],axis=1),3)
                A_aug = sym.Matrix(np.concatenate([temp_basis,vector],axis=1))
                red, _ = A_aug.rref(pivots=True)
                red = np.round(np.array(red.tolist(),dtype=np.float64),3)
                inconsistent = np.all(red[:,:-1] == 0,axis=1) & (red[:,-1] != 0)
                if not np.any(inconsistent):
                    res = red[:,-1]
                    integers = np.equal(np.mod(res, 1), 0)
                    lock_i = integers[0]
                    lock_j = integers[1]
                    lock_k = integers[2]
                if lock_i & lock_j & lock_k:
                    #print('A_aug\n',np.array(A_aug.tolist(),dtype=np.float64))
                    #print('Vector\n',vector)
                    #print('red\n',red)
                    #print('res\n',res)
                    #print('integers\n',integers)
                    return True
                else:
                    pass

    return False

def get_independent_basis(data,cutoff=100):
    independent_vectors = data.iloc[:1].values.reshape(1,3)
    basis = data.iloc[:1].values.reshape(1,3)
    cutoff = min(cutoff,len(data))
    for x in (range(1,cutoff,1)):
        temp_vector = data.iloc[x].values.reshape(1,3)
        # Check if any combination of basis leads to integer
        combination = check_basis_combinations(basis.T,temp_vector.T)
        if not combination:
            independent_vectors = np.append(independent_vectors,temp_vector,axis=0)
        basis = np.append(basis,temp_vector,axis=0)

    return independent_vectors


if __name__ == '__main__':
    main()