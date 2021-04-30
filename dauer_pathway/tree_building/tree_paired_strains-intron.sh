#!/bin/bash
#conda activate vcf-kit

for file in akt2.vcf daf2.vcf ;
do
bcftools view -Ov -S paired_strains.txt $file | bcftools view -Ov -e 'INFO/BCSQ ~ "intron"' \
> $file.paired_strains-intron.vcf
vk phylo tree nj $file.paired_strains-intron.vcf > $file.paired_strains-intron.newick
done
