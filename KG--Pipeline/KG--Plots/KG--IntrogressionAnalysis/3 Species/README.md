## **Introgression analysis using Dsuite (Malinsky et al. 2021)**

 This folder contains the Dsuite files and results used for D and F stats for the dataset considering four species (Kmar, Kher, and Ksp_ES).


### 1) Runnind Dsuite Dtrios

I ran Dsuite Dtrios to obtain D and F statistics for all possible species trios. As follows:
```
module load GCC/9.3.0
/scratch/waldirmbf/ES-Article_ABBABABA/Dstats/Dsuite/Build/Dsuite Dtrios /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_ReplacedNames.bcf /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/SETSfile_3species.txt -t /scratch/waldirmbf/ES-Article_ABBABABA/Dstats/KryptoTree_3species.newick -o 3species_Dtrios

```
On which:

`-t` = ['.newick' file](KryptoTree_3species.newick) with the species tree configuration supported by our phylogenetic reconstruction.

`-o` =  Prefix for output files. (i.e. 3species_Dtrios).

The server version used to run this command can be found at: [ES-Article-Dsuite_3species_DTRIOS.sbatch](ES-Article-Dsuite_3species_DTRIOS.sbatch).

Output files are:

[`3species_Dtrios_tree.txt`](3species_Dtrios_tree.txt) containing the trios according to the phylogenetic relantioships provided in the `-t` file.

[`3species_Dtrios_Dmin.txt`](3species_Dtrios_tree.txt) containing the trios according the minimum D for each trio, regardless of tree topology.

[`3species_Dtrios_BBAA.txt`](3species_Dtrios_BBAA.txt) containing the trios ordered by BBAA patterns (a potential attempt to infer populations/species relationships).

The [`combine.txt`](3species_Dtrios_combine.txt) and [`combine_stderr.txt`](3species_Dtrios_combine_stderr.txt) are output files used as inputs in `DTriosCombine` to look across specific  genomic regions.
