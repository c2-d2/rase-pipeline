# RASE prediction pipeline

## Introduction

This repository contains the RASE prediction pipeline. The method uses lineage
calling to identify antibiotic resistant clones from nanopore reads. In our
[paper](https://www.biorxiv.org/content/early/2018/08/29/403204), we
demonstrate on the example of pneumococcus that, using this approach,
antibiotic resistance can be predicted within minutes from the start of
sequencing.

## Overview of the pipeline

The main workflow is specified in `Snakefile`.  The pipeline starts by
detecting the provided RASE databases (in `database`) and MinION reads (in
`reads`). Both of these are then pre-processed: reads are sorted by time and the
database is uncompressed (i.e., the internal k-mer index reconstructed).

The reads are then compared to the database using
[ProPhyle](http://prophyle.github.io) and isolate, phylogroup, and resistance
to individual antibiotics predicted (as a function of time). Finally, the
obtained functions, as well as rank plots for selected moments, are plotted
using R.


## Installation of RASE

**Installing dependencies.**
Bioconda is the preferred way of installation of the software dependencies. We
recommend to create a separate software environment (here called `rase`):

```
conda create -n rase prophyle ete3 pysam snakemake samtools parallel r-optparse
```

and then activate it by

```
source activate rase
```

Alternatively, the packages can be installed into the default Conda environment.

```
conda install prophyle ete3 pysam snakemake samtools parallel r-optparse
```

Please note that, at some systems, the R package distributed by Conda might not
be properly built. The solution is then to create the environment without
`r-optparse`, and to install R and the Optparse package manually.


**Cloning the RASE pipeline.**
You can clone this repository using git

```
git clone https://github.com/c2-d2/rase-pipeline
```

or download it as a [single .tar.gz
file](https://github.com/c2-d2/rase-pipeline/archive/master.tar.gz).

**Installing a RASE database.** A RASE database should be placed into the
directory `database`.  Every database consists of two files: a compressed
ProPhyle index (`.tar`) and a table with metadata for individual database
isolates (`.tsv`). To be properly detected, both of the files should have the
same base name.

The default RASE database (Streptococcus pneumoniae, k=18) can be downloaded
from [RASE DB releases](https://github.com/c2-d2/rase-db/releases). A custom
database can be constructed using scripts and workflows from the [RASE DB
repository](https://github.com/c2-d2/rase-db).

**Placing nanopore reads.** Nanopore reads should be placed into the `reads`
directory as a single `.fq` file per sequencing experiment. Please, check the
suffix: `.fastq` files are not currently detected. Reads that were used in the paper
can be downloaded from https://zenodo.org/record/1405173.


## Running RASE

**Running prediction locally.** The RASE pipeline can be executed by the Make
command (all steps, i.e., preprocessing, predicting, and plotting):

```
make
```

This will run the entire RASE pipeline, including plotting. For the testing and
debugging purposes, it might be sometimes useful to run RASE on a small
example. The following command will detect the smallest provided database and
experiment, and will run RASE only for this single combination:

```
make test
```

**Running prediction on a cluster.** When multiple experiments and multiple
databases are combined, it can be useful to parallelize the computation. In the
default setting, RASE support the Harvard O2 cluster (Slurm-based), but the
configuration files are easy to adjust. Submitting to a cluster using Snakemake
can be done by `make o2`.

**Exporting outputs.** Outputs of the pipeline can be exported to a single
`.tar` archive by `make export`.

**Cleaning.** In some situations, it might be useful to clean intermediate
files.

* `make clean` - cleans all prediction outputs and plots, but keeps ProPhyle
  outputs (the most time consuming step)
* `make cleanall` - removes previous and also ProPhyle outputs
* `git clean -fxd` - removes all files that are not part of the git repository
  (including the databases and reads)


**Others.** For the list of available subcommands, see the output of `make
help`.


## Structure of the repository directories

* `benchmarks` - Snakemake benchmarks
* `database` - source database files
* `plots` - created plots
* `prediction` - intermediate prediction files
* `reads` - Minion reads
* `scripts` - RASE scripts
* `tests` - testing data


## Publication

Karel Brinda, Alanna Callendrello, Lauren Cowley, Themoula Charalampous, Robyn
S Lee, Derek R MacFadden, Gregory Kucherov, Justin O'Grady, Michael Baym,
William P Hanage. **Lineage calling can identify antibiotic resistant clones
within minutes.** [bioRxiv
403204](https://www.biorxiv.org/content/early/2018/08/29/403204), doi:
https://doi.org/10.1101/403204, 2018.

## Contact

Karel Brinda \<kbrinda@hsph.harvard.edu\>

