#!/bin/bash
#PBS -l select=1:ncpus=12:mem=48gb
#PBS -l walltime=00:30:00
##PBS -J 1-16

# Uncomment if running on HPC
#cd $PBS_O_WORKDIR

anaconda_path=~/anaconda3
env=CrystalAnalyser

source ${anaconda_path}/bin/activate
conda activate ${env}

mol2_directory=../PyCSD/molecular_packing_shells_20
n_splits=16
split_index=1
output_directory=geometric_interactions_packing_shells_20
output_prefix=''
n_processes=4

python geometric_interactions.py --mol2_directory=$mol2_directory --n_splits=$n_splits --split_index=$split_index --output_directory=$output_directory \
--output_prefix=$output_prefix --n_processes=$n_processes
