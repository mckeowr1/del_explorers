#!/bin/bash
#SBATCH -A b1059
#SBATCH -p b1059
#SBATCH -t 96:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="iqtree.oe"
#SBATCH --job-name="iqtree"

# Activate env
source activate	orthofinder_env

# Run MAFFT on [ORTHO_ID].fa
# Path to [ORTHO_ID].fa] can be full or relative

mafft --auto /projects/b1059/projects/Ryan/conservation/protein_fasta/OrthoFinder/Results_Jul20/Orthogroup_Sequences/OG0000001.fa > OG0000001.mafft

# Run IQTREE
# Path to [ORTHO_ID].mafft can be full or relative

iqtree -s OG0000001.mafft -nt 24

#This will generate a .treefile that you can quickly visualize in ITOL (https://itol.embl.de/upload.cgi)
