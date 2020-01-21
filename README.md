# Multi-tissue TWAS Simulation

Multi-tissue TWAS Simulation is an algorithm to simultaneously simulation genotypes, gene expression, and phenotype datasets across multiple tissues. Different from exisitng TWAS simnulation algorithms, the design of the Multi-tissue TWAS Simulation allows simulation of 1) polygenic genetic architecture of transcriptional regulation, 2) cross-tissue gene expression similarities, 3) flexible degrees of gene expression levels that are attributable to underlying regulatory eQTLs. 

## Prerequisites

The software is developed using R language and tested in Mac and Unix OS, and implemented in high performance computing (HPC) center. There is no required version of R, but R (>=3.4.2) is perferred. 

To run run_twas_simulation.R, the following libraries is required:
* optparse
* mvtnorm
* reshape2
* dplyr
* glmnet
* pls
* GBJ
* foreach
* doParallel

## Project Layout

The **software** folder contains the simulation and evaluation tools. 
* run_twas_simulation.R
  * The main script to generate simulated datasets
  * Automatically run several existing eQTL detecting methods and TWAS algorithms using the simulated datasets
* twas_simulation_util.R
  * Amenities of Multi-tissue TWAS Simulation, expected to be under the same folder as run_twas_simulation.R by default
* optim.so
  * Amenity of Multi-tissue TWAS Simulation, expected to be under the same folder as run_twas_simulation.R by default
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
* *--training 500* eQTL discovery dataset sample size (default = 500)
* *--testing 1000* TWAS sample size (default = 1000)
* *--snps 60* Number of SNPs in a given gene
* *--eqtls 30* Number of eQTLs among SNPs in a given gene 
* *--pct-mt-eqtls 0.8* Percentage of multi-tissue eQTLs among total eQTLs (rounded)
* *--maf 0.01,0.5* Range of minor allele frequency
* *--genes 100* Number of genes in a round of simulation
* *--expr-tis 10* Number of gene expressing tissues
* *--h2-ge 0.3* Heritability of gene expression levels
* *--cor-tissues 0.8* Similarities of gene expression levels among tissues, can be a set value like '0' for all tissue pairs, or a range like '-1,1'. In the latter case, Similarities of gene expression levels among tissues will be drawn from a uniform distribution ranging from -1 to 1. 
* *--r2-et 0.01* Variance of simulated traits explained by gene expression levels
* *--output-dir ./* Output directory of simulation results
* *--output-prefix twas_sim* Prefix of output files
* *--core 10* Number of cores to run parallel tasks
* *--random-seed 1* Random seed number

## Reference
The manuscript "Tissue specificity-aware TWAS framework identifies novel associations with metabolic and virologic traits in HIV-positive adults" is under preparation.

## Acknowledgements
We thank [Dr. Yiming Hu](https://github.com/Joker-Jerome/UTMOST), [Dr. Haky Im, and Alvaro N Barbeira](https://github.com/hakyimlab/MetaXcan) for their technical help and support.  
