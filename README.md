# Analysis workflow for "The genetic legacy of extreme exploitation in a polar vertebrate", Paijmans et al. Sci Rep 10, 5089 (2020)

## Overview

This is the workflow used in the paper ["The genetic legacy of extreme exploitation in a polar vertebrate"](https://doi.org/10.1038/s41598-020-61560-8).
The raw microsatellite data are available via the Zenodo repository: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4609701.svg)](https://doi.org/10.5281/zenodo.3585717)

To run the analyses, please download the complete folder. The folder named "Rcode" contains 14 scripts, which are named 1_ to 10_ . These scripts can be run a standard desktop machine.

Some of the analyses are computationally intensive. This is the case for the STRUCTURE analyses (scripts and data can be found in the folder "STRUCTURE") and the bottleneck analyses (scripts and data can be found in the folder "ABC"). The scripts in these folders should run on a server with sufficient computing power and memory.   

Most datasets which are produced along the way are already saved in subfolders, so that the analysis can be started at any point. In some cases the file size was too large and the file could not be included here. These files are available upon request.

In some cases, data needed to be transformed from STRUCTURE files to genepop files. This was done using [PGD Spider](http://www.cmpg.unibe.ch/software/PGDSpider/) and the .spid file is also given in this repository.

Finally, .bat files were used for mass conversion of files (see also script 9a_BOTTLENECK_in), also these .bat files are given in this repository.

## Basic statistics and PCA

### Within the folder "Rcode" the following scripts can be found:

- 1a_prep_data_ALL: prepping data for STRUCTURE
- 1b_hwe_ld_ALL: test for HWE and LD, removing loci that are not in HWE/LD
- 1f_pca_ALL: PCA (fig S12a & b)
- 2_remove_hybs: removing hybrids and A. tropicalis found using STRUCTURE
- 3a_prep_data_GAZ: prepping data for STRUCTURE, A. gazella only
- 4_stats: calculate # private alleles, Ar, Ho, Fis (input for table 1 and table S1)
- 5a_hwe: test for HWE (input for table S7)
- 5b_ld_test: test for LD
- 6_Fst_GAZ: calculate Fst (input for table S2)
- 7_pca_GAZ: PCA (fig S3a & b)
- 8_mratio_GAZ: calculate M ratio (input for table 1 and table S1) 
- 9a_BOTTLENECK_in: creates 1000 subsets of randomly drawn individuals for analysis with BOTTLENECK software
- 9b_BOTTLENECK_out: reads 1000 results of BOTTLENECK software and gets prop het ex values (input for fig 3 and fig S4)
- 10_sealing_effort: sealing effort (fig 1)

## STRUCTURE analysis

### Within the folder "STRUCTURE" the following scripts can be found:

(1) Folder "Trop_Gaz": in this STRUCTURE run we included all the data, so A. gazella, A. tropicalis and potential hybrids. We used a two population model (i.e. k = 2) to classify individuals that were admixed with at least 10% of the genetic attribution being to the secondary species (i.e. 0.10 ≤ q ≤ 0.9).
- 1c_run_structure: script to run STRUCTURE on the server
- 1d_parse_structure: parse and plot STRUCTURE output (fig S11)
- 1e_extract_hyb: get a list of hybrids identified as described above

(2) Folder "Trop_Gaz": in this STRUCTURE run we included only A. gazella individuals
- 3b_run_structure: script to run STRUCTURE on the server
- 3c_parse_structure: parse and plot STRUCTURE output (input fig 1, fig S1, S2)

## Coalescent-simulations and ABC analysis.
The scripts contained in the ABC folders were used to simulate genetic data under 
a bottleneck and a neutral demographic scenario using `strataG` as an interface
to [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/). These data
were then used for ABC analyses across genetic clusters/populations. 

### Prerequisites

(1) These script should be run on a multi-core machine, optimally using around 20
or more cores. However, everything can run quickly for testing purposes based
on a small number of simulations (say 1000 instead of 20000000)

(2) Install fastsimcoal2 and check using the [strataG](https://cran.r-project.org/web/packages/strataG/index.html) vignette
that the package can access fastsimcoal.

(3) Several other packages need to be installed, which are mentioned in the scripts.
Among them is a small package specifically written for the analysis of this paper,
the `sealABC` package, which can be install from GitHub with:

```
devtools::install_github("mastoffel/sealABC")
```

In addition, we slightly altered some specific functions of the package `sealABC`.
For this the script: mssumstatsAP is needed (in addition to the `sealABC` package). This script can be found in all folders where it is used.

### Within the folder "ABC" the following scripts can be found:

(1) Folder "fsc_cluster2019": simulations compared to emperical data on the level of genetic clusters
- a_prep_emp_data_5cluster: prepping emperical data
- b_sumstats_5cluster: calculating summary statistics from emperical data on genetic cluster level
- 1_simulate_diversity: simulates genetic diversity (nb the resulting file is not included on GitHub as the file size was too big. However it is available upon request)
- 2_ABCanalysis: ABC analysis part 1, Model selection and evaluation
- 3_ABCanalysis_posterior_distributions: ABC analysis part 2, Parameter estimation
-	4_abc_results: save ABC estimates to RData file (nb not all resulting files are included on GitHub as the file size was too big. However they are available upon request)
-	5_cv_eval_plots_FigS8: cross-validation plots (fig S8)
-	6_one_plot_Fig3: creates Fig 3 
-	7_posterior_predictive_checks: simulations for posterior predictive checks
-	8_posterior_predictive_plots_FigS7: creates figure for posterior predictive checks (Fig S7)
-	9_simulating_diversity_diff_Ne_Fig4a-c: simulations to calculate diversity loss using different Ne and plots (Fig 4a-c)
- 10_simulations_heatmap: simulations to calculate diversity loss using different Nebot size and duration
- 11_plotting_heatmap_Fig4d: creates heatmap (Fig 4d)
- Supp_conf_mat_plot_FigS6: confusion matrix plot (Fig S6)
-	Supp_density_plot_nehist_FigS9: density plots for Nehist (Fig S9)

(2) Folder "fsc_pop2019": simulations compared to emperical data on the level of populations (locality)
- a_prep_emp_data_8pop: prepping emperical data
- b_sumstats_8pop: calculating summary statistics from emperical data on genetic cluster level
- 1_simulate_diversity: simulates genetic diversity (nb the resulting file is not included on GitHub as the file size was too big. However it is available upon request)
- 2_ABCanalysis: ABC analysis part 1, Model selection and evaluation
- 3_ABCanalysis_posterior_distributions: ABC analysis part 2, Parameter estimation
-	4_abc_results: save ABC estimates to RData file (nb not all resulting files are included on GitHub as the file size was too big. However they are available upon request)
-	5_cv_eval_plots: cross-validation plots
-	6_one_plot_FigS4: creates Fig S4 

(3) Folder "sumstats_species": 
- sumstats_species: calculates summary statistics for other otariid species 
- compare_species_beeswarm_FigS5: creates plot showing allelic richness for all otariids (Fig S5)


The code is highly specific to the current analysis and probably has to be modified to be of use in other projects.
