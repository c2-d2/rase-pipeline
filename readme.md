# RASE prediction pipeline


## Introduction

This repository contains the RASE prediction pipeline. The method uses lineage
calling to identify antibiotic resistant clones from nanopore reads. In our
[paper](https://www.biorxiv.org/content/early/2018/08/29/403204), we
demonstrate on the example of pneumococcus that, using this approach,
antibiotic resistance can be predicted within minutes from the start of
sequencing. Please, look at the paper for more information.

The RASE [Snakemake](https://snakemake.readthedocs.io/) workflow is specified
within a single [Snakefile](Snakefile). When executed, the pipeline first
detects the provided RASE database(s) (in the `database` directory) and
nanopore reads (in the `reads` directory), and generates all
`<db>`-`<experiment>` combinations. In practice, the most common scenario is
usally "1 db vs. many experiments". After the detection step, reads and
database are pre-processed: reads are sorted by time of sequencing and the
database gets uncompressed (i.e., the full internal ProPhyle k-mer index
restored).  Subsequently, nanopore reads from individual experiments are
compared to the database(s) using [ProPhyle](http://prophyle.github.io), and
isolate, phylogroup and resistance to individual antibiotics predicted - all
as a function of time.  Finally, the obtained time characteristics, as well as
rank plots for selected moments, are visualized using R.


## Quick example

The following example demonstrates the power of the streptococcal RASE with metagenomic reads.
The entire computation takes only 6m on a standard laptop (MacBook Pro). Note
that this is the experiment from Figure 3 (with human reads removed in-silico).

To run the example, [install all
dependencies](https://github.com/c2-d2/rase/blob/master/environment.md) and run
the following code:


```bash
# clone this repository
git clone https://github.com/c2-d2/rase-predict

# download the default database
(cd rase-predict/database \
  && wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tar.gz \
  && wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tsv)

# download minion reads from a metagenomic experiment
(cd rase-predict/reads \
  && wget https://zenodo.org/record/1405173/files/sp10_norwich_P33.filtered.fq)

# run the pipeline
make -C rase-predict
```


## Installation

1) **Installing dependencies.** See [RASE computational
   enviroment](https://github.com/c2-d2/rase/blob/master/environment.md).

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

**Running prediction locally.** The RASE pipeline can be executed by the Make
command (all steps, i.e., preprocessing, predicting, and plotting):

```bash
make
```

This will run the entire RASE pipeline, including plotting. For the testing and
debugging purposes, it might be sometimes useful to run RASE on a small
example. The following command will detect the smallest provided database and
experiment, and will run RASE only for this single combination:

```bash
make test
```

**Running prediction on a cluster.** When multiple experiments and multiple
databases are combined, it can be useful to parallelize the computation. In the
default setting, RASE supports the [Harvard O2
cluster](https://rc.hms.harvard.edu/#cluster) (Slurm-based), but the
[configuration file](cluster.json) can be easily adjusted. Submitting to a
cluster using
[Snakemake](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution)
can be done by `make cluster`.

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


## Files and directories

* `benchmarks/` - Snakemake benchmarks. As soon as any step of the pipeline
  gets finished, a log file with information about timing and memory
  consumption will appear here.
* `database/` - Source database files.
   - Each database consists of two files: `<db>.tar.gz` and `<db>.tsv`.
* `logs/` - Logs from job submission systems.
* `plots/` - Plotted figures.
   - `<experiment>__<db>.timeline.pdf` - Prediction as a function of time.
   - `<experiment>__<db>.snapshot.<time>.pdf` - Rank plot for selected times (1
   minute, 5 minutes, last minute)
* `prediction/` - Prediction files.
   - `<experiment>__<db>.fq` - Nanopore reads after renaming and sorting by
   timestamp.
   - `<experiment>__<db>.bam` - Reads mapped onto the phylogenetic tree using
   ProPhyle.
   - `<experiment>__<db>/<timestamp>.tsv` - Cumulative weights calculated for
   all isolates from the DB at the time; `h1` corresponds to the weights used
   in the paper.
   - `<experiment>__<db>.predict.tsv` - Prediction timeline (each row
   corresponds to one minute).
* `reads` - Nanopore reads (`<experiment>.fq`).
* `scripts` - RASE scripts.
* `tests` - Testing data for scripts.


## FAQs

> Why am I getting 'libR.dylib Reason: image not found'?

On some systems, the R package distributed by Bioconda might not be properly
built and would display messages such as

```
dyld: Library not loaded: @rpath/libintl.9.dylib
   Referenced from: /Users/user/miniconda/envs/rase/lib/R/lib/libR.dylib
   Reason: image not found
Abort trap: 6
```

The solution is then to create the `raseenv` environment without `r-optparse`,
and to install R and the OptParse package manually.

> Why am I getting 'ETE: cannot connect to X server'?

ETE 3 library, which is used for tree plotting, internally depends on QT and
requires using an X-Server. This becomes problematic especially on virtual
machines.  For instance, on Ubuntu-based machines this can be solved by
installing several additional packages:

```
apt-get install xvfb libqt4-dev libgl1-mesa-dev libglu1-mesa-dev xauth xfonts-base
```

and then prepending the following string to commands for building the database.
```
xvfb-run --server-args="-screen 0 1024x768x24 -noreset" \
```


## Related repositories

* [RASE supplementary](http://github.com/c2-d2/rase-supplement). Supplementary Materials for the RASE paper, including figures and tables.
* [ProPhyle](http://prophyle.github.io). A highly accurate and resource-frugal DNA sequence classifier used by RASE.
* [Prophex](http://github.com/prophyle/prophex). A k-mer index based on the Burrows-Wheeler Transform, used by ProPhyle.


## License

[MIT](LICENSE).


## Contact

[Karel Brinda](https://scholar.harvard.edu/brinda) \<kbrinda@hsph.harvard.edu\>

