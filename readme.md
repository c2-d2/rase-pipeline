# RASE prediction pipeline

| This repository contains the RASE pipeline for rapid inference of antibiotic resistance and susceptibility using genomic neighbor typing using RASE. Other components of the software include the [RASE core package](https://github.com/c2-d2/rase/) and the [RASE DB skeleton](https://github.com/c2-d2/rase-db-skeleton). For more information, see the [RASE paper](https://www.biorxiv.org/content/10.1101/403204v2) and the [RASE supplementary materials](https://github.com/c2-d2/rase-supplement/). |
|-|


## Content

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Quick example](#quick-example)
* [Installation](#installation)
* [Running RASE](#running-rase)
  * [Subcommands](#subcommands)
* [Files and directories](#files-and-directories)
* [Related repositories](#related-repositories)
* [License](#license)
* [Contact](#contact)

<!-- vim-markdown-toc -->


## Introduction


This repository contains the RASE prediction pipeline for rapid inference of
antibiotic resistance and susceptibility using genomic neighbor typing. For
more information about the method and results, see the associated
[paper](https://www.biorxiv.org/content/10.1101/403204v2) and [supplementary
repository](https://github.com/c2-d2/rase-supplement).

The entire pipeline is specified within a single
[Snakemake](https://snakemake.readthedocs.io) [Snakefile](Snakefile), which can
be executed using GNU Make (see below).


## Quick example

The following example predicts antibiotic resistance from a sputum metagenomic
sample and the pneumococcal database.

After [install
dependencies](https://github.com/c2-d2/rase/blob/master/environment.md), run
the commands below. The entire computation should require
approximately 6 minutes on a standard laptop (MacBook Pro).


```bash
# clone and enter this repository
git clone --recursive https://github.com/c2-d2/rase-pipeline
# enter the directory
cd rase-pipeline
# run the pipeline on test data (S. pneumoniae database
# and nanopore reads from a metagenomic experiment)
make test
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

3) **Installing a RASE database.** A RASE database should be placed into
  `database/`.  Every database consists of two files: a compressed
  *k*-mer index (`.tar`) and a table with metadata for individual database
  strains (`.tsv`). To be properly detected, both of the files should have the
  same base name.

4) **Placing nanopore reads.** Nanopore reads should be placed into `reads/`
  as a single `.fq`/`.fastq`/`.fa`/`.fasta` file (possibly gzipped) per sequencing experiment.



## Running RASE

Upon execution using `make` the pipeline analyzes data in the following steps:

1) **Detection.** The pipeline detects user-provided RASE database(s) (in
`database/`) and read sets (in `reads/`), and generates all
`<db>`-`<reads>` experiments.

2) **Data preparation.** The pipeline decompressed the found databases and
sorts reads by time of sequencing. When the time information is not available
in the original reads (in read headers), it is estimated it based on the number
of processed basepairs (assuming constant flow (kbps per sec)).

3) **Matching.** Reads are matched against the databases using
[ProPhyle](https://prophyle.github.io/) and the computed nearest neighbors for
each read stored in `matching/` in the RASE-BAM format.

4) **Prediction.** Phenotypes are predicted from the computed nearest neighbors.

5) **Plotting.** The computed time characteristics of prediction are plotted.


### Subcommands

* `make` - Run everything.
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
* `reads/` - Nanopore reads (`<reads>.{fq,fastq,fa,fasta}{.gz,}`).

Output files:

* `matching/`
   - `<reads>.fa` - Nanopore reads sorted by time.
   - `<reads>__<db>.bam` - Read matches to the reference strains in RASE/BAM.
* `prediction/` - Prediction files.
   - `<reads>__<db>/<timestamp>.tsv` - Cumulative weights calculated for
   all strains from the database at the time; `h1` corresponds to the weights used in the paper.
   - `<reads>__<db>.predict.tsv` - Prediction timeline (each row
   corresponds to one minute).
* `plots/` - Plotted figures.
   - `<reads>__<db>.timeline.pdf` - Prediction as a function of time.
   - `<reads>__<db>.snapshot.<time>.pdf` - Rank plot for selected times (1
   minute, 5 minutes, last minute)


## Related repositories

* [*S. pneumoniae* RASE DB](https://github.com/c2-d2/rase-db-spneumoniae-sparc/).
* [*N. gonorrhoeae* RASE DB](http://github.com/c2-d2/rase-db-ngonorrhoeae-gisp).
* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

