# **Genetic analysis with SNPs from msGBS library**


### ES-Killfish Project Pipeline - by **George PACHECO** [![Foo](../KF_ES--GitHubAuxiliaryFiles/ORCIDGreenRoundIcon.png)](https://orcid.org/0000-0002-9367-6813) & **Waldir M. BERBEL-FILHO** [![Foo](../KF_ES--GitHubAuxiliaryFiles/ORCIDGreenRoundIcon.png)](https://orcid.org/0000-0001-6991-4685)


This documentation outlines the pipelines used for genetic analysis (SNPs genotypes and sites extracted from msGBS library) in the preprint manuscript

Please contact **Waldir Berbel-Filho** (waldirmbf@gmail.com) should any questions arise.
***
***

### 1) Data Access

 The raw msGBS files used for this manuscript data are stored at **SRA(XXXXX)** .

 The additional genome files incorporated in our analyses from (*Kryptolebias mamoratus* - DAN2K, HON9, LK1, and RHL; *Kryptolebias hermaphroditus* Central clade - Gitmo and PanRS) were extracted from [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188) (DAN2K, HON9, LK1, RHL and Gitmo) and [Soi et al. 2021](https://doi.org/10.1016/j.cbd.2020.100684) (PanRS). Following are the links to access these files:

| ID | Species | Population | Link | Reference |
| -------- | ---------- || --------  | -------- |
| DAN2K | *K. marmoratus* | Belize | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489803/SRR5489803 | [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188)
| HON9 | *K. marmoratus* | Honduras | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489795/SRR5489795 | [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188)
| LK1 | *K. marmoratus* | Florida | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489798/SRR5489798 | [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188)
| RHL | *K. marmoratus* | San Salvador | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489800/SRR5489800 | [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188)
| Gitmo | *K. hermaphroditus* - Central clade | Guantamo bay, Cuba | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489808/SRR5489808 | [Lins et al. 2018](https://doi.org/10.1139/gen-2017-0188)
| PanRS | *K. hermaphroditus* - Central clade| Panama | https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5489808/SRR5489808 | [Soi et al. 2021](https://doi.org/10.1016/j.cbd.2020.100684)

### 2) Genome Editing:

zcat ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.fna.gz | awk '{split($0,a," "); print a[1]'} > ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna
***

### 3) Genome Indexing (Bowtie2):

```
/data/home/waldir/Desktop/msGBS_data/Tools/bowtie2-2.3.4.3-linux-x86_64/bowtie2-2.3.5-linux-x86_64/bowtie2-build ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna ES-Genome
```
***

### 4) Genome Indexing (Samtools):

```
samtools faidx ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna
```
***

### 5) Mapping of msGBS samples:

```
while read i;
do
/data/home/waldir/Desktop/msGBS_data/Tools/bowtie2-2.3.4.3-linux-x86_64/bowtie2-2.3.5-linux-x86_64/bowtie2 -q --threads 8 -x ~/Desktop/msGBS_data/ES-Article/ES-Genome/ES-Genome -U ~/Desktop/msGBS_data/ES-Article/ES-GBS_Samples/$i.fq -S ~/Desktop/msGBS_data/ES-Article/ES-SAMs/$i.sam && samtools view -bS ~/Desktop/msGBS_data/ES-Article/ES-SAMs/$i.sam > ~/Desktop/msGBS_data/ES-Article/ES-First_BAMs/$i.bam && samtools sort ~/Desktop/msGBS_data/ES-Article/ES-First_BAMs/$i.bam -o ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam && samtools index ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam && mv ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam.bai ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bai
done < ~/Desktop/msGBS_data/ES-Article/ES-ListOfSamples_Edited_Strings.txt &> ~/Desktop/msGBS_data/ES-Article/ES-Mapping.txt
```
***


### 5) Reads processing and mapping for genome files | PaleoMix v1.3.2 #
```
module load Python
module load AdapterRemoval Bowtie2
module load SAMtools/1.10-GCC-8.3.0

python3 -m pip install paleomix==1.3.2 --user

paleomix bam dryrun --max-threads 6 --bowtie2-max-threads 3 --adapterremoval-max-threads 3 --log-file /scratch/waldirmbf/ES-Article_PaleoMix_Output/ES-Mapping-WGS_PaleoMix_Test.log --log-level info --destination /scratch/waldirmbf/ES-Article_PaleoMix_Output/ /scratch/waldirmbf/ES-Article_PaleoMix_Output/ES-Mapping-WGS_PaleoMix.yaml

```

###  6) Gets list of samples:

```
find /scratch/waldirmbf/ES-SortedIndexed/*.bam /scratch/waldirmbf/ES-Article_PaleoMix_Output/*.bam > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist
```
### 7) GLOBAL COVERAGE DISTRUBUTION | ANGSD--v0.929

```
/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((48*95/100)) -doCounts 1 -dumpCounts 2 -maxDepth $((48*1000)) -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.depth
```

##### Runs ANGSD for Sites:

```
/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((48*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((48*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES
```

##### _Number of SITES_: 115,397


##### Gets Real Coverage (_Genotype Likelihoods_):

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_SITES.GL-RealCoverage.txt
```

##### Gets Missing Data (_Genotype Likelihoods_):

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/115397}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_SITES.GL-MissingData.txt
```
#

***

### NJ PHYLOGENY | RAxML-NG--v1.0.1 ###

##### Converts the .haplo file into a .fasta file:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.haplo.gz | cut -f 4- | tail -n +2 | perl /home/waldirmbf/Software/Scripts/tsv_merge.pl --transp --ofs '' - | awk 'NR==FNR{id=$1; sub(".*\\/","",id); sub("\\..*","",id); x[FNR]=id} NR!=FNR{ print ">"x[FNR]"\n"$1}' /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.labels - > /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES.fasta

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_SITES.sbatch

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_NoHybrids_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_NoHybrids_SITES.sbatch

# Here we use RAxML-NG to search for a ML Phylogeny based on this FASTA alingment: 10 Random PLUS 10 Parsimony:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --search --tree pars{100},rand{100} --model GTR+G --site-repeats on --log PROGRESS --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES.fasta --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_SITES.sbatch


# Then, we use RAxML-NG to bootstrap the BEST ML Phylogeny:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --bootstrap --model GTR+G --bs-trees 100 --site-repeats on --log PROGRESS --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES.fasta --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.raxml.bestTree --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.BOOTs

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_SITES.sbatch

# Finally, we add the bootstrap values supports to the generated ML phylogeny:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 2 --support --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.raxml.bestTree --bs-trees /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.BOOTs.raxml.bootstraps --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.SUPPORTED

***

###                         ###
# SNP CALLING | ANGSD--v0.929 #
###                         ###

>>>> Dataset III (AllSamples -- NoKgraNoKbra (33) / SITES):

# BAM List:

grep -v "Kgra\|Kbra" /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist

# ANGSD Run:

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist  -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES

> # of SITES: 863,662

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36943160_stderr.txt

#### Gets Real Coverage:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_NoKbraNoKgra_SITES.GL-RealCoverage.txt

##### Gets Missing Data:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/863662}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_NoKbraNoKgra_SITES.GL-MissingData.txt

>>>> Dataset III (AllSamples -- NoKgraNoKbra (33) / SNPs):

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.04 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs

> # of SNPs: 9,532

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SNPs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SNPs.sbatch

>>> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36942918_stderr.txt

##### Gets Real Coverage:

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.GL-RealCoverage.txt
```

##### Gets Missing Data:

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/9532}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.GL-MissingData.txt
```


### ESTIMATION OF INDIVIDUAL ANCESTRIES | ngsAdmix

```
export N_REP=100

for K in `seq -w 2 3`
do
    NGSADMIX_BIN=/home/waldirmbf/Software/ngsAdmix/NGSadmix/NGSadmix /home/waldirmbf/Software/Scripts/wrapper_ngsAdmix.sh -P 7 -likes /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz -K $K -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o /scratch/waldirmbf/ES-Article_ngsAdmix/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.${K}

done

chmod +x /scratch/waldirmbf/ES-Article_ngsAdmix/ToRunadmix_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ngsAdmix/ToRunadmix_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch

cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_ngsAdmix/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.annot
```
***

###  MULTIDIMENSIONAL SCALING | ngsDist + get_PCA.R

## Here are perform a multidimensional scaling anlyse on the genetic distance matrix created above:

##### Gets genetic distance matrix:

```
install.packages('libgsl')

/home/waldirmbf/Software/ngsDist/ngsDist --n_threads 2 --geno /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.beagle.gz --pairwise_del --seed 11 --probs --n_ind 48 --n_sites 1810 --labels /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels --out /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.dist

/home/waldirmbf/Software/ngsDist/ngsDist --n_threads 2 --geno /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz --pairwise_del --seed 18 --probs --n_ind 33 --n_sites 9532 --labels /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels --out /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.dist
```

##### Performs MDS:

```
tail -n +3 /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.dist | Rscript --vanilla --slave /home/waldirmbf/Software/Scripts/get_PCA.R --no_header --data_symm -n 10 -m "mds" -o /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.mds

tail -n +3 /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.dist | Rscript --vanilla --slave /home/waldirmbf/Software/Scripts/get_PCA.R --no_header --data_symm -n 10 -m "mds" -o /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.mds
```

##### Creates `.annot` file:

```
cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.annot

cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.annot
```
***

### HETEROZYGOSITY CALCULATION

## Here we calculate the percentage of heterozygous genotypes in our NoSNPCalling sites.

##### Generates a `.bed` file based on the `.mafs` file:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.bed

##### Creates a `.pos` file based on this new `.bed` and index it accordingly:

```
awk '{print $1"\t"($2+1)"\t"$3}' /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.bed > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.pos


/home/waldirmbf/Software/angsd/angsd sites index /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.pos
```

##### Gets files:

```
parallel --plus --dryrun /home/waldirmbf/Software/angsd/angsd -i {} -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -anc /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -sites /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.pos -GL 1 -doSaf 1 -fold 1 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/{/...} :::: /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES.sbatch

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_GBSs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_GBSs.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/jobname_37857280_stderr.txt

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ToRunHet_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/jobname_37869116_stderr.txt
> /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/jobname_37870855_stderr.txt
```

##### Gets fractions:

```
parallel --plus "/home/waldirmbf/Software/angsd/misc/realSFS {} > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/{/..}.het" ::: /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/*.saf.idx
```

##### Calculates the percentage of heterozygous sites:

```
fgrep '.' *.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | gawk '{print $1"\t"$2"\t"lol[1]}' | sort -k 1,1gr | awk '{split($0,a,"."); print a[1]"\t"$2"\t"$3'} > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.Heterozygosity.txt
```

##### These results were plotted using the Rscript below:

PBGP--ToPlotProportionOfHeterozygousSites.R
***
#
#
#

### Gets VCF

/home/waldirmbf/Software/ANGSD-30/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_ABBABABA/ES-Article--AllSamples_NoKbraNoKgra_SITES_TEST.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((3*95/100)) -GL 1 -doPost 1 -doMajorMinor 1 -doMaf 1 -doVcf 1 -doGeno 1 -doCounts 1 -doGlf 2 -MinMaf 0.04 -SNP_pval 1e-6 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((3*600)) -dumpCounts 2 -postCutoff 0.95 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/Test/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_TEST

### Gets Geno-Depth File

tail -n +5253 ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.bcf | less -S
