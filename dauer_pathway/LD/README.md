So there are 3 scripts.

The 1st cegwas_genotype_matrix.sh takes a list of alleles (Regions file format for Bcftools) and spits out a genotype martrix where ref and alt values have been replaced by 1 and -1. I have it set up as a loop rn to go through multiple lists of allele pairs.

The 2nd LD_paired_alleles.R takes the genotype matrices and spits out a matrix with R^2 values between the alleles. Some of it I adapted from Nemascan so I’m pretty confident I know what it does but happy to help you dig into it more. So those two will get you the pairwise LD between variants.

The 3rd script plot_LD I have is just to make a heatmap R^2 values, the LD_paired_alleles.R just spits out half the matrix so I wanted to get rid of the NA in the top half.
