

# **Genetic analysis with SNPs from msGBS library**
ES-Killfish Project Pipeline | Ultra Documentation - by George Pacheco [![Foo](../ORCID-iD.png)](https://orcid.org/0000-0002-9367-6813)  and Waldir M. Berbel-Filho [![Foo](../ORCID-iD.png)](https://orcid.org/0000-0001-6991-4685)


This documentation outlines the pipelines used for genetic analysis (SNPs genotypes and sites extracted from msGBS library) in the preprint manuscript

Last Modified :XXXXXX            

Please, contact george.pacheco@snm.ku.dk or waldirmbf@gmail.com should any question arise.

__________________________________________

####               ###
# KFGP -- ES-Article #
###   4 Oct 2020   ###

# Genome Editing:

zcat ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.fna.gz | awk '{split($0,a," "); print a[1]'} > ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna

# Genome Indexing (Bowtie2):

/data/home/waldir/Desktop/msGBS_data/Tools/bowtie2-2.3.4.3-linux-x86_64/bowtie2-2.3.5-linux-x86_64/bowtie2-build ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna ES-Genome

# Genome Indexing (Samtools):

samtools faidx ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna

# Mapping:

while read i;
do
/data/home/waldir/Desktop/msGBS_data/Tools/bowtie2-2.3.4.3-linux-x86_64/bowtie2-2.3.5-linux-x86_64/bowtie2 -q --threads 8 -x ~/Desktop/msGBS_data/ES-Article/ES-Genome/ES-Genome -U ~/Desktop/msGBS_data/ES-Article/ES-GBS_Samples/$i.fq -S ~/Desktop/msGBS_data/ES-Article/ES-SAMs/$i.sam && samtools view -bS ~/Desktop/msGBS_data/ES-Article/ES-SAMs/$i.sam > ~/Desktop/msGBS_data/ES-Article/ES-First_BAMs/$i.bam && samtools sort ~/Desktop/msGBS_data/ES-Article/ES-First_BAMs/$i.bam -o ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam && samtools index ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam && mv ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bam.bai ~/Desktop/msGBS_data/ES-Article/ES-SortedIndexed/$i.bai
done < ~/Desktop/msGBS_data/ES-Article/ES-ListOfSamples_Edited_Strings.txt &> ~/Desktop/msGBS_data/ES-Article/ES-Mapping.txt

###                                          ###
# GLOBAL COVERAGE DISTRUBUTION | ANGSD--v0.929 #
###                                          ###

~/Desktop/msGBS_data/Tools/ngsTools/angsd/angsd -nThreads 2 -ref ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fna -bam ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((42*95/100)) -doCounts 1 -dumpCounts 2 -maxDepth $((42*1000)) -out ~/Desktop/msGBS_data/ES-Article/ES-ANGSDRuns/ES-Article--AllSamples.depth

> # of SITES: 257,984

###                         ###
# SNP CALLING | ANGSD--v0.929 #
###                         ###

~/Desktop/msGBS_data/Tools/ngsTools/angsd/angsd -nThreads 2 -ref ~/Desktop/msGBS_data/ES-Article/ES-Genome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((42*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.03 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((42*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -doVcf 1 -out ~/Desktop/msGBS_data/ES-Article/ES-ANGSDRuns/ES-Article--AllSamples_SNPs

> # of SNPs: 2,060

# Real Coverage Calculation:

zcat ~/Desktop/msGBS_data/ES-Article/ES-ANGSDRuns/ES-Article--AllSamples_SNPs.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.labels - > ~/Desktop/msGBS_data/ES-Article/ES-Miscellaneous/ES-RealCoverage/ES-Article--AllSamples_SNPs.GL-RealCoverage.txt

# Missing Data Calculation:

zcat ~/Desktop/msGBS_data/ES-Article/ES-ANGSDRuns/ES-Article--AllSamples_SNPs.beagle.gz | tail -n +2 | perl ~/Desktop/msGBS_data/Tools/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.labels - | awk '{print $1"\t"$3"\t"$3*100/2060}' > ~/Desktop/msGBS_data/ES-Article/ES-Miscellaneous/ES-MissingData/ES-Article--AllSamples_SNPs.GL-MissingData.txt

###                                              ###
#  MULTIDIMENSIONAL SCALING | ngsDist + get_PCA.R  #
###                                              ###

## Here are perform a multidimensional scaling anlyse on the genetic distance matrix created above:

# To get distance matrix:

~/Desktop/msGBS_data/Tools/ngsTools/ngsDist/ngsDist --n_threads 8 --geno ~/Desktop/msGBS_data/ES-Article/ES-ANGSDRuns/ES-Article--AllSamples_SNPs.beagle.gz --pairwise_del --seed 33 --probs --n_ind 42 --n_sites 2060 --labels ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.labels --out ~/Desktop/msGBS_data/ES-Article/ES-MDS/ES-Article--AllSamples_SNPs.dist

# To perform MDS:

tail -n +3 ~/Desktop/msGBS_data/ES-Article/ES-MDS/ES-Article--AllSamples_SNPs.dist | Rscript --vanilla --slave ~/Desktop/msGBS_data/Tools/Scripts/get_PCA.R --no_header --data_symm -n 10 -m "mds" -o ~/Desktop/msGBS_data/ES-Article/ES-MDS/ES-Article--AllSamples_SNPs.mds

# Create .annot file:

cat ~/Desktop/msGBS_data/ES-Article/ES-Lists/ES-Article--AllSamples.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > ~/Desktop/msGBS_data/ES-Article/ES-MDS/ES-Article--AllSamples_SNPs.annot

###                                            ###
# ESTIMATION OF INDIVIDUAL ANCESTRIES | ngsAdmix #
###                                            ###

export N_REP=100

for K in `seq -w 2 10`
do
    echo /groups/hologenomics/fgvieira/scripts/wrapper_ngsAdmix.sh -P 4 -likes ~/data/Temp/Files/ES-Article--AllSamples_SNPs.beagle.gz -K $K -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o ~/data/Temp/ES-Admix/ES-Article--AllSamples_SNPs.${K}

done | xsbatch -c 4 --mem-per-cpu 1024 --max-array-jobs 9 -J ngsAdmix -R --time 10-00 --

###         ###
# PHYLOGENIES #
###         ###

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR5489795.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR5489798.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR5489800.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR5489803.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR5489808.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

> ~/Desktop/msGBS_data/Tools/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/SRR9139883.1 --split-files --skip-technical -F --outdir ~/Desktop/msGBS_data/ES-Article/ES-OtherGenomes/FASTQs/ --gzip

###                                             ###
# READS' PROCESSING AND MAPPING | PaleoMix v1.3.2 #
###                                             ###

module load Python
module load AdapterRemoval Bowtie2
module load SAMtools/1.10-GCC-8.3.0

python3 -m pip install paleomix==1.3.2 --user

paleomix bam dryrun --max-threads 6 --bowtie2-max-threads 3 --adapterremoval-max-threads 3 --log-file /scratch/waldirmbf/ES-Article_PaleoMix_Output/ES-Mapping-WGS_PaleoMix_Test.log --log-level info --destination /scratch/waldirmbf/ES-Article_PaleoMix_Output/ /scratch/waldirmbf/ES-Article_PaleoMix_Output/ES-Mapping-WGS_PaleoMix.yaml

paleomix bam run --max-threads 6 --bowtie2-max-threads 6 --adapterremoval-max-threads 6 --log-file /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/ES-Mapping-WGS_OutGroup_PaleoMix.log --log-level info --destination /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/ /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/ES-Mapping-WGS_OutGroup_PaleoMix.yaml

> less /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/jobname_37497312_stderr.txt

chmod +x /scratch/waldirmbf/ES-Article_PaleoMix_Output/ToRunPaleoMix.sbatch
chmod +x /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/ToRunPaleoMix_WGS_OutGroup.sbatch

sbatch /scratch/waldirmbf/ES-Article_PaleoMix_Output/ToRunPaleoMix.sbatch
sbatch /scratch/waldirmbf/ES-Article_PaleoMix_OutGroup_Output/ToRunPaleoMix_WGS_OutGroup.sbatch

###                                          ###
# GLOBAL COVERAGE DISTRUBUTION | ANGSD--v0.929 #
###                                          ###

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((48*95/100)) -doCounts 1 -dumpCounts 2 -maxDepth $((48*1000)) -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.depth

> # of SITES: 256,247

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.depth.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.depth.sbatch

###                                           ###
# CREATION OF SPECIFIC DATASETS | ANGSD--v0.932 #
###                                           ###

>>> Dataset 0 | ALL GOOD SAMPLES (SITES / 48 SAMPLES):

# BAM List:

find /scratch/waldirmbf/ES-SortedIndexed/*.bam /scratch/waldirmbf/ES-Article_PaleoMix_Output/*.bam > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist

# ANGSD Run:

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((48*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((48*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES

> # of SITES: 115,397

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36812980_stderr.txt

# Real Coverage Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_SITES.GL-RealCoverage.txt

# Missing Data Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/115397}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_SITES.GL-MissingData.txt

>>> Dataset 1 | ALL GOOD SAMPLES / NO Hybrids (SITES / 46 SAMPLES):

# BAM List:

fgrep -v -f /scratch/waldirmbf/ES-Article_Lists/ES-Article--BadSamples.list /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.BAMlist

# ANGSD Run:

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.BAMlist  -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((46*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((46*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES

> # of SITES: 115,909

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36904983_stderr.txt

# Real Coverage Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_NoHybrids_SITES.GL-RealCoverage.txt

# Missing Data Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/115909}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_NoHybrids_SITES.GL-MissingData.txt

### NJ PHYLOGENY | RAxML-NG--v1.0.1 ###

# First we convert the HAPLO file into a FASTA:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoHybrids_SITES.haplo.gz | cut -f 4- | tail -n +2 | perl /home/waldirmbf/Software/Scripts/tsv_merge.pl --transp --ofs '' - | awk 'NR==FNR{id=$1; sub(".*\\/","",id); sub("\\..*","",id); x[FNR]=id} NR!=FNR{ print ">"x[FNR]"\n"$1}' /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoHybrids_SITES.labels - > /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES.fasta

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_SITES.sbatch

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_NoHybrids_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToConvertHAPLOintoFASTA_ES-Article--AllSamples_NoHybrids_SITES.sbatch

# Here we use RAxML-NG to search for a ML Phylogeny based on this FASTA alingment: 10 Random PLUS 10 Parsimony:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --search --tree pars{100},rand{100} --model GTR+G --site-repeats on --log PROGRESS --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES.fasta --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_SITES.sbatch

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --search --tree pars{100},rand{100} --model GTR+G --site-repeats on --log INFO --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES.fasta --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_NoHybrids_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToRunUnguidedMLPhylogeny_ES-Article--AllSamples_NoHybrids_SITES.sbatch

# Then, we use RAxML-NG to bootstrap the BEST ML Phylogeny:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --bootstrap --model GTR+G --bs-trees 100 --site-repeats on --log PROGRESS --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES.fasta --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.raxml.bestTree --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.BOOTs

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_SITES.sbatch

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 10 --bootstrap --model GTR+G --bs-trees 100 --site-repeats on --log INFO --msa /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES.fasta --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s.raxml.bestTree --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s.BOOTs

> chmod +x /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_NoHybrids_SITES.sbatch
> sbatch /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ToBootBESTPhy_ES-Article--AllSamples_NoHybrids_SITES.sbatch

# Finally, we add the bootstrap values supports to the generated ML phylogeny:

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 2 --support --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.raxml.bestTree --bs-trees /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.BOOTs.raxml.bootstraps --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_SITES_100s.SUPPORTED

/home/waldirmbf/Software/RAxML-NG/raxml-ng --threads 2 --support --tree /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s.raxml.bestTree --bs-trees /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s.BOOTs.raxml.bootstraps --prefix /scratch/waldirmbf/ES-Article_Phylogenies/ES-Article_ML/ES-Article--AllSamples_NoHybrids_SITES_100s.SUPPORTED

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###                         ###
# SNP CALLING | ANGSD--v0.929 #
###                         ###

>>>> Dataset II (AllSamples):

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((48*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.03 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((48*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs

> # of SNPs: 1,810

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.sbatch

# Real Coverage Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_WithWGSs_SNPs.GL-RealCoverage.txt

# Missing Data Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/1810}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_WithWGSs_SNPs.GL-MissingData.txt

>>>> Dataset III (AllSamples -- NoKgraNoKbra (33) / SITES):

# BAM List:

grep -v "Kgra\|Kbra" /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist

# ANGSD Run:

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist  -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES

> # of SITES: 863,662

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.sbatch

> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36943160_stderr.txt

# Real Coverage Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_NoKbraNoKgra_SITES.GL-RealCoverage.txt

# Missing Data Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/863662}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_NoKbraNoKgra_SITES.GL-MissingData.txt

>>>> Dataset III (AllSamples -- NoKgraNoKbra (33) / SNPs):

/home/waldirmbf/Software/angsd/angsd -nThreads 2 -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -bam /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.BAMlist -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -minInd $((33*95/100)) -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -MinMaf 0.04 -SNP_pval 1e-6 -doPost 2 -doGeno 3 -doPlink 2 -geno_minDepth 3 -setMaxDepth $((33*600)) -dumpCounts 2 -postCutoff 0.95 -doHaploCall 1 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs

> # of SNPs: 9,532

chmod +x /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SNPs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SNPs.sbatch

>>> /scratch/waldirmbf/ES-Article_ANGSDRuns/jobname_36942918_stderr.txt

# Real Coverage Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article_ES-RealCoverage/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.GL-RealCoverage.txt

# Missing Data Calculation:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/9532}' > /scratch/waldirmbf/ES-Article_Miscellaneous/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.GL-MissingData.txt

###                                            ###
# ESTIMATION OF INDIVIDUAL ANCESTRIES | ngsAdmix #
###                                            ###

export N_REP=100

for K in `seq -w 2 3`
do
    NGSADMIX_BIN=/home/waldirmbf/Software/ngsAdmix/NGSadmix/NGSadmix /home/waldirmbf/Software/Scripts/wrapper_ngsAdmix.sh -P 7 -likes /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz -K $K -minMaf 0 -tol 1e-6 -tolLike50 1e-3 -maxiter 10000 -o /scratch/waldirmbf/ES-Article_ngsAdmix/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.${K}

done

chmod +x /scratch/waldirmbf/ES-Article_ngsAdmix/ToRunadmix_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch
sbatch /scratch/waldirmbf/ES-Article_ngsAdmix/ToRunadmix_ES-Article--AllSamples_WithWGSs_NoKbraNkgra_SITES_WGSs.sbatch

cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_ngsAdmix/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.annot

###                                              ###
#  MULTIDIMENSIONAL SCALING | ngsDist + get_PCA.R  #
###                                              ###

## Here are perform a multidimensional scaling anlyse on the genetic distance matrix created above:

# To get distance matrix:

install.packages('libgsl')

/home/waldirmbf/Software/ngsDist/ngsDist --n_threads 2 --geno /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_WithWGSs_SNPs.beagle.gz --pairwise_del --seed 11 --probs --n_ind 48 --n_sites 1810 --labels /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels --out /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.dist

/home/waldirmbf/Software/ngsDist/ngsDist --n_threads 2 --geno /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.beagle.gz --pairwise_del --seed 18 --probs --n_ind 33 --n_sites 9532 --labels /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels --out /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.dist

# To perform MDS:

tail -n +3 /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.dist | Rscript --vanilla --slave /home/waldirmbf/Software/Scripts/get_PCA.R --no_header --data_symm -n 10 -m "mds" -o /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.mds

tail -n +3 /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.dist | Rscript --vanilla --slave /home/waldirmbf/Software/Scripts/get_PCA.R --no_header --data_symm -n 10 -m "mds" -o /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.mds

# Create .annot file:

cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_WithWGSs_SNPs.annot

cat /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_NoKbraNoKgra_SITES.labels | awk '{split($0,a,"_"); print $1"\t"a[1]"_"a[3]}' > /scratch/waldirmbf/ES-Article_MDS/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs.annot

#### HETEROZYGOSITY CALCULATION ####

## Here we calculate the percentage of heterozygous genotypes in our NoSNPCalling sites.

# First we generate a '.bed' file based on the '.mafs' of this run:

zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_NoKbraNoKgra_SITES.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.bed

# After we create a position file based on this new  '.bed' and index it accordingly usings ANGSD:

awk '{print $1"\t"($2+1)"\t"$3}' /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.bed > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.pos


/home/waldirmbf/Software/angsd/angsd sites index /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.pos

# Getting files:

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

# Getting fractions:

parallel --plus "/home/waldirmbf/Software/angsd/misc/realSFS {} > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/{/..}.het" ::: /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/*.saf.idx

# Finally, we calculate the percentage of heterozygous sites:

fgrep '.' *.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | gawk '{print $1"\t"$2"\t"lol[1]}' | sort -k 1,1gr | awk '{split($0,a,"."); print a[1]"\t"$2"\t"$3'} > /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article_Het/ES-Article--AllSamples_NoKbraNoKgra_SITES.Heterozygosity.txt

# We locally plot these results using the Rscript below:

> PBGP--ToPlotProportionOfHeterozygousSites.R
