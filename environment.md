# RASE computational environment

The easiest way how to setup a computational environment for RASE is using
[Bioconda](https://bioconda.github.io/). This approach has been tested on
multiple Unix and OS X machines, including clusters and virtual machines.

## Dependencies

* [Python 3](https://www.python.org/downloads/)
* [ProPhyle](http://prophyle.github.io)
* [ETE 3](http://etetoolkit.org/)
* [PySAM](https://github.com/pysam-developers/pysam)
* [GNU Make](https://www.gnu.org/software/make/) or equivalent
* [GNU parallel](https://www.gnu.org/software/parallel/)
* [Ghost Script](https://www.ghostscript.com/)
* [Pandas](https://pandas.pydata.org/)
* [SnakeMake](https://snakemake.readthedocs.io)
* [SAMtools](http://www.htslib.org/)
* [R](https://www.r-project.org/)
* [R OptParse](https://cran.r-project.org/web/packages/optparse/)
* [GCC 4.8+](https://gcc.gnu.org/) or equivalent
* [zlib](https://zlib.net/)


## Setting up an environment

### Using a separate BioConda environment

We recommend to create a separate software environment (here called `rase`):

```bash
conda create -n rase \
	prophyle ete3 pysam snakemake-minimal samtools parallel r-optparse pandas
```

The environment can then be activated by

```bash
source activate rase
```

### Using the default BioConda environment

Alternatively, the packages can also be installed directly into the default
BioConda environment. Nevertheless, this is not always reliable since some of
the RASE dependencies might collide with packages that were installed
previously.

```bash
conda install prophyle ete3 pysam snakemake samtools parallel r-optparse
```

### Alternative ways of installation

All the dependencies can also be installed without BioConda.

Many of these packages are distributed using standard package systems such as
[APT](https://wiki.debian.org/Apt).

```bash
apt-get install build-essential python3 zlib1g-dev r-base r-cran-optparse ghostscript
```

All the Python packages (ProPhyle, PySAM, ETE 3, and Snakemake) can be
installed using [PIP](https://pypi.org/project/pip/):

```bash
pip3 install prophyle pysam ete3 snakemake
```


## Known issues

### libR.dylib Reason: image not found

At some systems, the R package distributed by Bioconda might not be properly
built and would display messages such as

```
dyld: Library not loaded: @rpath/libintl.9.dylib
  Referenced from: /Users/user/miniconda/envs/rase/lib/R/lib/libR.dylib
  Reason: image not found
Abort trap: 6
```

The solution is then to create the `rase` environment without `r-optparse`, and
to install R and the OptParse package manually.

### ETE: cannot connect to X server

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

