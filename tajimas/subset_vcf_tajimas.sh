#!/bin/bash
#SBATCH -J h1a      ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 02:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=16G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1            ## Number of processors for your task


#bcftools view -Oz /projects/b1059/analysis/WI-20210121/isotype_only/WI.20210121.hard-filter.isotype.bcsq.20210401.vcf.gz \
#-r III -f TYPE="snp" > chromIII.vcf


for regions in III III:2994514-III:3040846 X X:14147700-X:14152070 IV IV:420011-IV:425177 \
I I:10750498-I:10776703 ;
do
bcftools view -Ou /projects/b1059/analysis/WI-20210121/isotype_only/WI.20210121.hard-filter.isotype.bcsq.20210401.vcf.gz \
-r $regions -f TYPE="snp" |
bcftools annotate -Oz -x FORMAT,INFO > $regions.vcf.gz #Must remove annotations so that we can run through LD python script
done
