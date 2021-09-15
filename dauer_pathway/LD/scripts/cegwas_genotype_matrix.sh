#!/bin/bash
#SBATCH -J geno_matrix_loop      ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 12:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=24G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1            ## Number of processors for your task


#Script to go from VCF to Genotype Matrix for DAF-2 AKT-2 Paired alleles
#CHROM POS REF ALT

#module load gatk
module load bcftools

for filename in /projects/b1059/projects/Ryan/dauer/paired_alleles/*_alleles ;

do



bcftools view \
-R $filename \
-Ou /projects/b1059/projects/Ryan/WI.20210121.hard-filter.isotype.vcf.gz |
#> paired_allele.vcf

#gatk VariantsToTable \
#-V paired_allele.vcf \
#-F CHROM -F POS -F REF -F ALT -GF GT  \
#-O paired_allele_gm.tsv

bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
  sed 's/[[# 0-9]*]//g' | \
  sed 's/:GT//g' | \
  sed 's/0|0/-1/g' | \
  sed 's/1|1/1/g' | \
  sed 's/0|1/NA/g' | \
  sed 's/1|0/NA/g' | \
  sed 's/.|./NA/g'  | \
  sed 's/0\/0/-1/g' | \
  sed 's/1\/1/1/g'  | \
  sed 's/0\/1/NA/g' | \
  sed 's/1\/0/NA/g' | \
  sed 's/.\/./NA/g' > $filename.Genotype_Matrix.tsv


done
