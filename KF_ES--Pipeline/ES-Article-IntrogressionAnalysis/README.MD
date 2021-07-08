## **Introgression analysis using Dsuite (Malinsky & Svardal, 2021)**


This folder contains all files used in the Introgression analysis ran using the software **Dsuite v--0.4 r28** [Malinsky et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13265).

### 1) Generating .bcf file as input to Dsuite

I had to reran the SNP calling pipeline in **ANGSD--v0.9.3**, to create a .bcf file as input for the Dsuite analysis. For that, I used the following command:  
```
/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_ABBABABA/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.04 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doBcf 1 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs
```
The server version used to run this command can be found at: [ES-Article-SNPcalling_doVcF.sbatch](ES-Article-SNPcalling_doVcF.sbatch).

### 1) Generating .bcf file as input to Dsuite
I manually edited the names in the original .bcf file to the new names used in [ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_ReplacedNames.bcf](ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_ReplacedNames.bcf).

This file was then used as the input VCF file for the Dsuite analysis.
