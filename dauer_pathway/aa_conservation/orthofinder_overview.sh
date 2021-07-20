#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="orthofinder.oe"
#SBATCH --job-name="orthofinder"


### build your env
# $ conda create -n orthofinder_env
# $ source activate orthofinder_env
# $ conda install -c bioconda orthofinder
# $ conda install -c bioconda mafft
# $ conda install -c bioconda iqtree

# Activate env
source activate orthofinder_env

# Run Orthofinder
# Path to [FASTA_DIR] can be full or relative
# [FASTA_DIR] contains all your sample FASTAs + GOI FASTAs

orthofinder -f [FASTA_DIR] -og -t 12

# Orthofinder results will be in ~/[FASTA_DIR]/OrthoFinder/Results[DATE]
# ~/Results[DATE]/Orthogroups/Orthogorups.txt contains all orthogroup IDs and the gene IDs within each orthogroup
# ~/Results[DATE]/Orthogroup_Sequences/ contains all orthogroup FASTAs
# Parse through Orthogroups.txt by eye, and identify the orthogroup ID [ORTHO_ID] that contains your gene of interest


# You will need to run the next steps in a separate script (unless you want to automatize the previous manual step!)



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
