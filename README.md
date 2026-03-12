<p align="right"> <img src="phylothin.jpg" width="200"> </p>

# *PhyloThin* 

This github repository provides the software *PhyloThin*.

With *PhyloThin*, we present a fully automated statistical tool that can detect and correct for sampling bias in prokaryotic populations. Removing the strong sampling bias of bacteria allows to identify the effective amount of information and prevents misleading biased conclusions in various analysis, such as the estimation of the pangenome size and gene frequencies.

## Step 1: Install R

- [R](https://www.R-project.org/)

```
sudo apt install r-base
```

### Required Packages for R

- ape
- phangorn
- dplyr

```
R
install.packages("ape")
install.packages("phangorn")
install.packages("dplyr")
q()
```

### References

1. R Core Team 2021. R: A language and environment for statistical computing. R Foundation for
  Statistical Computing, Vienna, Austria. https://www.R-project.org/
2. Paradis E. and Schliep K. 2019. ape 5.0: An Environment for Modern Phylogenetics and Evolutionary
Analyses in R. Bioinformatics 35: 526-528.
3. Schliep K.P. 2011. phangorn: Phylogenetic Analysis in R. Bioinformatics 27(4): 592-593.
4. Hadley Wickham, Romain François, Lionel Henry, Kirill Müller and Davis Vaughan 2023. dplyr: A
  Grammar of Data Manipulation. https://CRAN.R-project.org/package=dplyr

## Step 2: Clone GitHub-Repository

```
git clone https://github.com/fbaumdicker/phylothin.git
cd phylothin
```

## Step 3: Additional Requirements

- [PATHd8](https://www2.math.su.se/PATHd8/)

*PhyloThin* uses *PATHd8* as default for transforming the input tree to an ultrametric tree (all tips are equidistant from the root). For this *PATHd8* needs to be saved as *PATHd8* in the same folder as the input files (see input *path_to_folder*). [Here](https://www2.math.su.se/PATHd8/) you can find already compiled versions of *PATHd8*. Otherwise compile *PATHd8* on your own as follows.

```
cd pathd8_source
cc PATHd8.c -O3 -lm -o PATHd8
mv PATHd8 ../PATHd8
```

If you prefer your own tool (or the tree is already ultrametric), *PATHd8* is not needed and can be skipped (i.e., no need to install *PATHd8*) by adding *"no_PATHd8"* at the end of the command-line input to run *PhyloThin*.

### Reference

1. Tom Britton, Cajsa Lisa Anderson, David Jacquet, Samuel Lundqvist, and Kåre Bremer 2007.
Estimating Divergence Times in Large Phylogenetic Trees. Systematic Biology 56(5): 741–752.

## Step 4: Run *PhyloThin*

```
Rscript phylothin.r path_to_folder input_tree (priority_list) (no_PATHd8) (no_clade)
```

### Input

- *path_to_folder*: path to the folder where the input is stored and the output-folder will be generated.
- *input_tree*: file-name of the inferred phylogenetic tree (rooted, not necessary ultrametric) stored as a Newick-file in *path_to_folder*. As a first step the input tree is made ultrametric by using *PATHd8*. If you prefer your own tool to make the tree ultrametric, do so and give the ultrametric tree as the input tree.
- *priority_list*: optional csv-file stored in *path_to_folder* and containing the sample ids (same as the tip labels of the input tree) in the first column and their priority in the second column. Preferably, samples with priority *0* are removed, then those without priority *NA*, then those with low priority (high number) until finally those with the highest priority *1* are removed.
- *no_PATHd8*: If *input_tree* is already ultrametric (you may have used your preferred tool) and *PATHd8* is not installed, you can skip *PATHd8* in *PhyloThin* by adding *"no_PATHd8"* in the command-line input to run *PhyloThin*.\
!! CAUTION !! Make sure that the given *input_tree* is already ultrametric !! Otherwise *PhyloThin* may give you wrong results.
- *no_clade*: If no clade-output is prefered (see output clades.csv below) add *"no_clade"* in the command-line input to run *PhyloThin*. This may save some running time of *PhyloThin*. For trees with more than 10,000 tips, the clade-ouput is not possible.

### Output

The input tree name is always used as a suffix in all output file names.

Folder *phylothinoutput/*: every output of *PhyloThin* is stored in this folder.
- *clades.csv*: table with sample ids and the oversampled clade to which they belong. The first column gives the sample id, the second column indicates if the sample is part of a oversampled clade (same number for samples in the same clade, the number itself has no further meaning). *NA* in the clade column means that the sample is not part of an oversampled clade. Therefore, each clade consists of removed samples and one (sometimes more are possible) sample which has been kept by *PhyloThin*. For further analysis it may sometimes also be meaningful to keep all samples and weight them according to how many samples belong to the same oversampled clade. For trees with more than 10,000 tips, this clade-ouput is not possible.
- *removed_ids.txt*: list of sample ids which should be removed for further analysis.
- *kept_ids.txt*: list of sample ids which should be kept for further analysis.
- *reduced_tree.nwk*: the (reduced) phylogenetic tree where biased samples are removed.
- *treecomparison.pdf*: comparison of the full and reduced phylogenetic trees. Removed samples are marked with red crosses in the full tree.
- *tree.pdf*: plot of the full tree where tips are labeled and removed samples are marked with red crosses.

Folder *phylothinoutput/check/*: additional visualizations which may help to understand the decisions of *PhyloThin*. See *Supplementary Note 1* in the accompanying manuscript for further details.
- *treecutting.pdf*
- *scaling.pdf*
- *um_treecomparison.pdf*: comparison of the full and reduced ultrametric, phylogenetic trees. Removed samples are marked with red crosses in the full ultrametric tree.

### Large Phylogenetic Trees - Subsampling

For very large sample sizes the mutation rate may not be high enough to allow a reliable inference of the shortest coalescent times in the phylogenetic tree. To prevent issues with *PhyloThin*, we suggest to subsample the strain collection multiple times. By considering the ratio of subsamples where a strain has been identified as oversampled we can then identify oversampling in the complete data set (see the accompanying manuscript for further details). This can be done as follows:

```
Rscript phylothin_subsampling.r path_to_folder input_tree (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)
```

For details on the computations of the default parameter setting, see "Supplementary Note 1: PhyloThin – Implementation" in the accompanying manuscript.

### Low mutation rate - Mutation-Sensitive Thinning

For low mutation rates, branches of length zero become more and more common (especially in large phylogenetic trees) which may violate the assumptions of *PhyloThin* and may lead to an overestimation of sampling bias. We therefore recommend ensuring that a sufficient number of variable sites were used to infer the phylogenetic tree. If this is not possible, we suggest cross-checking the results obtained with *PhyloThin* using the mutation-sensitive version of *PhyloThin*. In this version, an inferred mutation rate is incorporated into the statistical test.
*Mutation-Sensitive Thinning* can be performed by providing the additional input *number_variable_sites* (the number of SNP positions used to infer the phylogenetic tree), as follows:

```
Rscript phylothin.r path_to_folder input_tree (priority_list) (no_PATHd8) (no_clade) -m number_variable_sites
```

CAUTION: *Mutation-Sensitive Thinning* may underestimate sampling bias if the resolution of the phylogenetic tree is insufficient. To mitigate this effect, *Mutation-Sensitive Thinning* can be combined with *Subsampling* as follows:

```
Rscript phylothin_subsampling.r path_to_folder input_tree (-r number_of_subsamples) (-s subsample_size) (no_PATHd8) -m number_variable_sites
```

## Test - Examples

```
Rscript phylothin.r ./Example ListeriaMonocytogenes.nwk
```

```
Rscript phylothin.r ./Example um_ListeriaMonocytogenes.nwk no_PATHd8
```


**Date: 12.03.2026**
