#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem=16000
#SBATCH --job-name=ANGSD_SNPs_doVCF
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=20:00:00
#SBATCH --mail-user=waldirmbf@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/scratch/waldirmbf/ES-Article_ABBABABA
################################################

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_ABBABABA/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.04 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doBcf 1 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs
