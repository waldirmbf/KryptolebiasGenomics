#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem=16000
#SBATCH --job-name=Dsuite_4species_Dtrios
#SBATCH --output=jobname_%J_stdout.txt
#SBATCH --error=jobname_%J_stderr.txt
#SBATCH --time=20:00:00
#SBATCH --mail-user=waldirmbf@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/scratch/waldirmbf/ES-Article_ABBABABA/Dstats
################################################
module load GCC/9.3.0
/scratch/waldirmbf/ES-Article_ABBABABA/Dstats/Dsuite/Build/Dsuite Dtrios /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_ReplacedNames.bcf /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/SETSfile.txt -t /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/KryptoTree.newick -o 4species_Dtrios
