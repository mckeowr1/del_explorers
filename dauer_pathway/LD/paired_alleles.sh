#Dont use this script, use  Genotype Matrix from Cegwas



#Script to go from VCF to Genotype Matrix for DAF-2 AKT-2 Paired alleles
#CHROM POS REF ALT


module load gatk
module load bcftools

bcftools view \
-R paired_alleles.tsv \
-Ov /projects/b1059/projects/Ryan/WI.20210121.hard-filter.isotype.vcf.gz \
-e GT="mis"
> paired_allele.vcf

gatk VariantsToTable \
-V paired_allele.vcf \
-F CHROM -F POS -F REF -F ALT -GF GT  \
-O paired_allele_gt.tsv
