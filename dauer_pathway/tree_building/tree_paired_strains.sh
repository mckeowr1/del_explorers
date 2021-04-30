#!/bin/bash
#conda activate vcf-kit

for file in akt2.vcf daf2.vcf ;
do
bcftools view -Ov -S paired_strains.txt $file > $file.paired_strains.vcf
vk phylo tree nj $file.paired_strains.vcf > $file.paired_strains.newick
done
