#!/bin/bash

#SBATCH --job-name=UK_0420_1221
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5GB
#SBATCH --qos=standard
#SBATCH --time=5-00:00:00

python "scripts/Cross_neutralization/Compute_FR.py" "UK_0420_1221/results/SpikeGroups.pck" "UK_0420_1221/results/Mutation_Profiles.pck" "UK_0420_1221/results/epitope_data/dms_per_ab_per_site.csv" "ALL" "None" "UK_0420_1221/Cross_with_delta.pck" "UK_0420_1221/results/Cross_react_dic_spikegroups_ALL.pck" "cluster_True" 50
