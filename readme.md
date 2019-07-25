# RASE prediction pipeline


## Introduction

This repository contains the RASE prediction pipeline. The workflow is
specified within a single [Snakemake](https://snakemake.readthedocs.io)
[Snakefile](Snakefile) can be executed using GNU Make (see below).


## Quick example

The following example predicts antibiotic resistance from a sputum metagenomic
sample and the pneumococcal database.

After [install
dependencies](https://github.com/c2-d2/rase/blob/master/environment.md), run
the commands below. The entire computation should require
approximately 6m on a standard laptop (MacBook Pro).


```bash
# clone and enter this repository
git clone --recursive https://github.com/c2-d2/rase-pipeline
cd rase-pipeline

# download the default database
wget -P database https://github.com/c2-d2/rase-db-spneumoniae-sparc/releases/download/v1.3/spneumoniae-sparc.k18.{tsv,tar.gz}

# download minion reads from a metagenomic experiment
wget -P reads wget https://zenodo.org/record/1405173/files/sp10_norwich_P33.filtered.fq

# run the pipeline
make
```


## Installation

1) **Installing dependencies.** See [RASE computational
   enviroment](environment.md).

2) **Cloning the RASE prediction pipeline.**
    You can clone this repository using git

    ```bash
    git clone --recursive https://github.com/c2-d2/rase-pipeline
    ```

  or download it as a [single .tar.gz
  file](https://github.com/c2-d2/rase-predict/archive/master.tar.gz).

3) **Installing a RASE database.** A RASE database should be placed into the
  directory `database`.  Every database consists of two files: a compressed
  ProPhyle index (`.tar`) and a table with metadata for individual database
  isolates (`.tsv`). To be properly detected, both of the files should have the
  same base name.

  The default RASE database (Streptococcus pneumoniae, k=18) can be downloaded
  from [RASE DB releases](https://github.com/c2-d2/rase-db/releases). A custom
  database can be constructed using scripts and workflows from the [RASE DB
  repository](https://github.com/c2-d2/rase-db).

4) **Placing nanopore reads.** Nanopore reads should be placed into the `reads`
  directory as a single `.fq` file per sequencing experiment. Please, check the
  suffix: `.fastq` files are not currently detected. Also, the pipeline assumes
  that the provided reads keep the original naming convention from ONT. Reads
  that were used in the paper can be downloaded from
  https://zenodo.org/record/1405173.


## Running RASE

Upon execution using `make` the pipeline analyzes data in the following steps:

*1) Detection.* The pipeline detects user-provided RASE database(s) (in
`database/`) and read sets (in `reads/`), and generates all
`<db>`-`<reads>` experiments.

*2) Data preparation.* The pipeline decompressed the found databases and sorts
reads by time of sequencing. When the time information is not available in the
original reads, it is estimated it based on the number of processed basepairs
(assuming constant flow (kbps per sec)).

*3) Matching.* Reads are matched against the databases using
[ProPhyle](https://prophyle.github.io/) and the computed nearest neighbors for
each read stored in `matching/` in the RASE-BAM format.

*4) Prediction.* Phenotypes are predicted from the computed nearest neighbors.

*5) Plotting.* The computed time characteristics of prediction are plotted.


### Subcommands

* `make`, `make all` - Run everything.
* `make cluster` - Submit jobs to a cluster.
* `make export` - Export all outputs to `rase_results.tar`.
* `make clean` - Clean plots and prediction files.
* `make cleanall` - Clean all output files (including bam files and logs).
* `make replot` - Re-plot all figures.
* `make test` - Run the smallest experiment only. If database or reads are not present, download examples. For testing and debugging purposes.


## Files and directories

Input files:
* `database/` - Source database files.
   - Each database consists of two files: `<db>.tar.gz` and `<db>.tsv`.
* `reads` - Nanopore reads (`<reads>.fq`).

Output files:

* `matching/`
   - `<reads>.fa` - Nanopore reads sorted by time.
   - `<reads>__<db>.bam` - Read matches to the reference strains in RASE/BAM.
* `prediction/` - Prediction files.
   - `<reads>__<db>/<timestamp>.tsv` - Cumulative weights calculated for
   all isolates from the DB at the time; `h1` corresponds to the weights used
   in the paper.
   - `<reads>__<db>.predict.tsv` - Prediction timeline (each row
   corresponds to one minute).
* `plots/` - Plotted figures.
   - `<reads>__<db>.timeline.pdf` - Prediction as a function of time.
   - `<reads>__<db>.snapshot.<time>.pdf` - Rank plot for selected times (1
   minute, 5 minutes, last minute)


## Related repositories

* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

