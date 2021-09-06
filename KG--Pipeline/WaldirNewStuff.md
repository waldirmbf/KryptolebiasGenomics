

##### Gets a `.labels` file:

```
awk '{split($0,a,"/"); print a[5]}' /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist | awk '{split($0,b,"."); print b[1]}' > /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels
```

##### Gets Real Coverage (_Genotype Likelihoods_):

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.counts.gz | tail -n +2 | gawk ' {for (i=1;i<=NF;i++){a[i]+=$i;++count[i]}} END{ for(i=1;i<=NF;i++){print a[i]/count[i]}}' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - > /scratch/waldirmbf/ES-Article_ANGSDRuns/RealCoverage/ES-Article--AllSamples_SITES.GL-RealCoverage.txt
```

##### Gets Missing Data (_Genotype Likelihoods_):

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.beagle.gz | tail -n +2 | perl /home/waldirmbf/Software/Scripts/call_geno.pl --skip 3 | cut -f 4- | awk '{ for(i=1;i<=NF; i++){ if($i==-1)x[i]++} } END{ for(i=1;i<=NF; i++) print i"\t"x[i] }' | paste /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.labels - | awk '{print $1"\t"$3"\t"$3*100/115397}' > /scratch/waldirmbf/ES-Article_ANGSDRuns/MissingData/ES-Article--AllSamples_SITES.GL-MissingData.txt
```
#



##### Generates a `.bed` file based on the `.mafs` file:

```
zcat /scratch/waldirmbf/ES-Article_ANGSDRuns/ES-Article--AllSamples_SITES.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - > /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_SITES.bed
```

##### Creates a position file based on this new `.bed`:

```
awk '{print $1"\t"($2+1)"\t"$3}' /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_SITES.bed > /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_SITES.pos
```

##### Indexs the `.pos` file created above:

```
/home/waldirmbf/Software/angsd/angsd sites index /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_SITES.pos
```

##### Gets files:

```
parallel --plus --dryrun /home/waldirmbf/Software/ANGSD-30/angsd/angsd -i {} -ref /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -anc /home/waldirmbf/ES-Article_REFGenome/GCA_007896545.1_ASM789654v1_genomic.Edited.fasta -sites /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_SITES.pos -GL 1 -doSaf 1 -fold 1 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 30 -minQ 20 -out /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/{/...} :::: /scratch/waldirmbf/ES-Article_Lists/ES-Article--AllSamples_SITES.BAMlist > /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ToRunHet_ES-Article--AllSamples_NoHybrids_SITES.sbatch
```

##### Gets fractions:

```
parallel --plus "/home/waldirmbf/Software/ANGSD-30/angsd/misc/realSFS {} > /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/{/..}.het" ::: /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/*.saf.idx
```

##### Calculates the percentage of heterozygous sites:

```
fgrep '.' *.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | gawk '{print $1"\t"$2"\t"lol[1]}' | sort -k 1,1gr | awk '{split($0,a,"."); print a[1]"\t"$2"\t"$3'} > /scratch/waldirmbf/ES-Article_ANGSDRuns/Het/ES-Article--AllSamples_NoHybrids_SITES.Heterozygosity.txt
```