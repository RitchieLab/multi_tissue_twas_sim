# Multi-tissue TWAS Simulation

Multi-tissue TWAS Simulation is an algorithm to simultaneously simulation genotypes, gene expression, and phenotype datasets across multiple tissues. Different from exisitng TWAS simnulation algorithms, the design of the Multi-tissue TWAS Simulation allows simulation of 1) polygenic genetic architecture of transcriptional regulation, 2) cross-tissue gene expression similarities, 3) flexible degrees of gene expression levels that are attributable to underlying regulatory eQTLs. 

## Prerequisites

The software is developed using R language and tested in Mac and Unix OS, and implemented in high performance computing (HPC) center. There is no required version of R, but R (>=3.5.0) is perferred. 

To run simulation using run_twas_simulation.R, the following libraries are required:
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
* [mvtnorm](https://cran.r-project.org/web/packages/mvtnorm/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
* [dplyr](https://dplyr.tidyverse.org/articles/dplyr.html)
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

The following packages are needed for eQTL detection and TWAS 
* [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
* [pls](https://cran.r-project.org/web/packages/pls/index.html)
* [GBJ](https://cran.r-project.org/web/packages/GBJ/index.html)


## Project Layout

The **software** folder contains the simulation and evaluation tools. 
* run_twas_simulation.R
  * The main script to generate simulated datasets
  * Automatically run several existing eQTL detecting methods and TWAS algorithms using the simulated datasets
* twas_simulation_util.R
  * Amenities of Multi-tissue TWAS Simulation, expected to be under the same folder as run_twas_simulation.R by default
* optim.*
  * Amenity of eQTL detection using the group LASSO algorithm. Expected to be under the same folder as run_twas_simulation.R by default
* run_h2_est.R
  * An evaluation tool to estimate genetic heritability of simulated traits
  
## Setup and Example using Command Lines

The following example assumes that you have installed all prerequisite.

1) Clone the repository.
```bash
$ git clone https://github.com/RitchieLab/multi_tissue_twas_sim
```

2) Go to the software folder.
```bash
$ cd multi_tissue_twas_sim/software
```

3) Run the Simulation script, run_twas_simulation.R.
```bash
$ Rscript run_twas_simulation.R \
  --training 500 \
  --testing 1000 \
  --snps 60 \
  --eqtls 30 \
  --pct-mt-eqtls 0.8 \
  --maf 0.01,0.5 \
  --genes 100 \
  --expr-tis 10 \
  --h2-ge 0.3 \
  --cor-tissues 0.8 \
  --r2-et 0.01 \
  --output-dir ./ \
  --output-prefix twas_sim \
  --core 10 \
  --random-seed 1
```

The example command parameters include all that one will need for multi-tissue TWAS simulation. Each specific parameter means the followings:
* *--training* eQTL discovery dataset sample size (default = 500)
* *--testing* TWAS sample size (default = 1000)
* *--snps* Number of SNPs in a given gene
* *--eqtls* Number of eQTLs among SNPs in a given gene 
* *--pct-mt-eqtls* Percentage of multi-tissue eQTLs among total eQTLs (rounded)
* *--maf* Range of minor allele frequency
* *--genes* Number of genes in a round of simulation
* *--expr-tis* Number of gene expressing tissues
* *--h2-ge* Heritability of gene expression levels
* *--cor-tissues* Similarities of gene expression levels among tissues, can be a set value like '0' for all tissue pairs, or a range like '-1,1'. In the latter case, Similarities of gene expression levels among tissues will be drawn from a uniform distribution ranging from -1 to 1. 
* *--r2-et* Variance of simulated traits explained by gene expression levels
* *--output-dir* Output directory of simulation results
* *--output-prefix* Prefix of output files
* *--core* Number of cores to run parallel tasks
* *--random-seed* Random seed number
* *--simulation* Only generate simulated genotype, gene expression, and phenotype datasets. Do not run eQTL detection and TWAS.

***Please note that the above command line also runs only the simulation part of the algorithm. To run subsequent eQTL detecting methods and TWAS algorithms that are embeded in the script, remove the flag "--simulation".***

## Embedded eQTL Detecting methods
1) [Elastic Net](https://www.nature.com/articles/ng.3367) as implemented in PrediXcan
2) [Group LASSO](https://www.nature.com/articles/s41588-019-0345-7) as implemented in UTMOST

## Embedded TWAS Algorithms
1) single-tissue gene-level association, e.g. PrediXcan ([link to GitHub](https://github.com/hakyimlab/MetaXcan), [link to paper](https://www.nature.com/articles/ng.3367))
2) cross-tissue gene-level association, e.g. MulTiXcan ([link to GitHub](https://github.com/hakyimlab/MetaXcan), [link to paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007889)), UTMOST ([link to GitHub](https://github.com/Joker-Jerome/UTMOST), [link to paper](https://www.nature.com/articles/s41588-019-0345-7))

## Reference
The manuscript "Tissue specificity-aware TWAS framework identifies novel associations with metabolic and virologic traits in HIV-positive adults" is under review.

## Acknowledgements
We thank [Dr. Yiming Hu](https://github.com/Joker-Jerome/UTMOST), [Dr. Haky Im, and Alvaro N Barbeira](https://github.com/hakyimlab/MetaXcan) for their technical help and support.  
