#!/bin/bash
#PBS -l select=1:ncpus=12:mem=48gb
#PBS -l walltime=12:00:00
##PBS -J 1-16

# Uncomment if running on HPC
#cd $PBS_O_WORKDIR

anaconda_path=~/anaconda3
env=CrystalAnalyser

source ${anaconda_path}/bin/activate
conda activate ${env}

input_directory=ex_geometric_interaction_dir
n_shells=12
output_directory=neighbour_shell_descriptors_test
output_prefix=''
n_processes=10

python get_neighbour_shell_descriptors.py --input_directory=$input_directory --n_shells=$n_shells --output_directory=$output_directory \
--output_prefix=$output_prefix --n_processes=$n_processes
