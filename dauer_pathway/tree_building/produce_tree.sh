#!/bin/bash
#conda activate vcf-kit

for file in akt2.vcf daf2.vcf ;
do
vk phylo tree nj $file > $file.newick
done
