# RASE prediction pipeline

## Introduction

This repository contains the RASE prediction pipeline. The method uses lineage
calling to identify antibiotic resistant clusters from MinION reads. In our
paper, demonstrate on the example of pneumococcus that, using this approach,
antibiotic resistance can be predicted within minutes.

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

The default RASE database (S.pneumoniae, k=18) can be downloaded from [RASE DB
releases](https://github.com/c2-d2/rase-db/releases). A custom database can be
constructed using scripts and workflows from the [RASE DB
repository](https://github.com/c2-d2/rase-db).

**Placing nanopore reads.** Nanopore reads should be placed into the `reads`
directory as one `.fq` for each experiment.


## Running RASE

**Running prediction locally.** The RASE pipeline can be executed by the Make command:
```
make
```

This will run the entire RASE pipeline, including plotting.

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

## Contact

Karel Brinda \<kbrinda@hsph.harvard.edu\>

