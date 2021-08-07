## **Introgression analysis using Dsuite (Malinsky et al., 2021)**


Given the genetic structure found between samples of the selfing *Kryptolebias* from Central and South America (most likely to represent different lineages within the same species, see [Lira et al. 2021](https://onlinelibrary.wiley.com/doi/epdf/10.1111/jfb.14753)), on this dataset we treat the samples from Panama and Southwest Cuba as a different population/clade ('Kher_CENTRAL'), as strongly supported by our phylogenetic analysis.

Therefore, this folder contains the Dsuite files and results used for D and F stats for the dataset considering four species (Kmar, Kher_CENTRAL, Kher_SOUTH, and Ksp_ES).


### 1) Runnind Dsuite Dtrios

I ran Dsuite Dtrios to obtain D and F statistics for all possible species trios. As follows:
```
module load GCC/9.3.0
/scratch/waldirmbf/ES-Article_ABBABABA/Dstats/Dsuite/Build/Dsuite Dtrios /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_ReplacedNames.bcf /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/SETSfile.txt -t /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/KryptoTree.newick -o 4species_Dtrios
```
On which:

`-t` = ['.newick' file](KryptoTree.newick) with the species tree configuration supported by our phylogenetic reconstruction.

`-o` =  Prefix for output files. (i.e. 4species_Dtrios).

The server version used to run this command can be found at: [ES-Article-Dsuite_4species_DTRIOS.sbatch](ES-Article-Dsuite_4species_DTRIOS.sbatch).

Output files are:

[`4species_Dtrios_tree.txt`](4species_Dtrios_tree.txt) containing the trios according to the phylogenetic relantioships provided in the `-t` file.

[`4species_Dtrios_Dmin.txt`](4species_Dtrios_tree.txt) containing the trios according the minimum D for each trio, regardless of tree topology.

[`4species_Dtrios_BBAA.txt`](4species_Dtrios_BBAA.txt) containing the trios ordered by BBAA patterns (a potential attempt to infer populations/species relationships).

The [`combine.txt`](4species_Dtrios_combine.txt) and [`combine_stderr.txt`](4species_Dtrios_combine_stderr.txt) are output files used as inputs in `DTriosCombine` to look across specific  genomic regions.
