#!/bin/bash
#conda activate vcf-kit

for file in akt2.vcf daf2.vcf ;
do
bcftools view -Ov -e 'INFO/BCSQ ~ "intron"' $file > $file.-intron.vcf
vk phylo tree nj $file.-intron.vcf > $file.-intron.newick
done
