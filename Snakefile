"""
Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import datetime
import glob
import itertools
import os
import re
import sys

mj, mn = sys.version_info.major, sys.version_info.minor
if mj < 3 or mn < 6:
    print("!!!!", file=sys.stderr)
    print("!!!! RASE PIPELINE ERROR  ", file=sys.stderr)
    print("!!!!", file=sys.stderr)
    print("!!!! Python 3.6+ is required ({}.{} provided)".format(mj, mn), file=sys.stderr)
    print("!!!!", file=sys.stderr)
    sys.exit(1)

#
# 0) Config
#

snakemake.shell.prefix("set -eo pipefail;")
configfile: "config.yaml"
localrules: all
#localrules: plot, plot_timeline, plot_snapshots, preprocess_reads, all

re_read_files = re.compile(r'(.*)/([^/].*?)\.(fasta|fastq|fa|fq)(\.gz)?', re.IGNORECASE)
re_index = re.compile(r'(.*)/([^/].*?)\.tar\.gz', re.IGNORECASE)

#
# 1) Detect indexes
#


def file_size(fn):
    return os.stat(fn).st_size


def smallest_file(fns):
    size = 10**30
    sfn = None
    for fn in fns:
        rfn = os.path.realpath(fn)
        if file_size(rfn) < size:
            size = file_size(rfn)
            sfn = fn
    return sfn


print(file=sys.stderr)
print("RASE DB and reads detection is starting", file=sys.stderr)
print(file=sys.stderr)

indexes = [x for x in glob.glob("database/*.tar.gz")]
if len(indexes) == 0:
    print("!!!! ", file=sys.stderr)
    print("!!!! RASE PIPELINE ERROR ", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    print("!!!! No index provided (db.{tsv,tar.gz}).", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    databases = []
    smallest_index = None
    smallest_database = None
else:
    databases = [re_index.match(x).group(2) for x in indexes]
    smallest_index = sorted(indexes, key=file_size, reverse=True)[0]
    smallest_database = re_index.match(smallest_index).group(2)
    print("  RASE databases:    ", ", ".join(databases), file=sys.stderr)
    print("  Smallest database: ", smallest_database, file=sys.stderr)
    print(file=sys.stderr)

#
# 2) Detect reads
#
readfiles = [x for x in glob.glob("reads/*") if re_read_files.search(x)]
if len(readfiles) == 0:
    print("!!!! ", file=sys.stderr)
    print("!!!! RASE PIPELINE ERROR ", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    print("!!!! No files with sequencing reads provided", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    experiments = []
    smallest_readfile = None
    smallest_experiment = None
else:
    experiments = [re_read_files.match(x).group(2) for x in readfiles]
    smallest_readfile = sorted(readfiles, key=file_size, reverse=True)[0]
    smallest_experiment = re_read_files.match(smallest_readfile).group(2)
    print("  Reads:             ", ", ".join(experiments), file=sys.stderr)
    print("  Smallest dataset:  ", smallest_experiment, file=sys.stderr)
    print(file=sys.stderr)

#
# 3) Check dependencies
#

snakemake.shell("(../scripts/rase_test_environments.sh 2>&1) >/dev/null || ../scripts/rase_test_environments.sh")

rule all:
    input:
        [
            [
                [
                    ancient(f"matching/.{e}.read.complete"),
                    ancient(f"matching/.{e}__{db}.assign.complete"),
                    ancient(f"prediction/.{e}__{db}.predict.complete"),
                    ancient(f"plots/.{e}__{db}.plot.complete"),
                ]
                for db in databases
            ]
            for e in experiments
        ],

rule database:
	input:
		[ancient(f"database/.{db}.complete") for db in databases]


rule test:
    input:
        [
            [
                [
                    ancient(f"matching/.{e}.read.complete"),
                    ancient(f"prediction/.{e}__{db}.predict.complete"),
                    ancient(f"plots/.{e}__{db}.plot.complete"),
                ]
                for db in [smallest_database]
            ]
            for e in [smallest_experiment]
        ],


for suffix in [f"{x}{y}" for x in ("fa", "fq", "fasta", "fastq") for y in ("", ".gz")]:
    rule:
        input:
            reads=ancient("reads/{pref}."+suffix)
        output:
            t="matching/.{pref}.read.complete"
        params:
            reads="matching/{pref}.fa"
        benchmark:
            "benchmarks/{pref}.readprep.log"
        shell:
            """
                minion=1 && ../src/rase/rase_has_datetime.py "{input.reads}" || minion=0; \
                echo "Has datetimes: $minion";
                if [ "$minion" -eq "1" ]; then
                    ../src/rase/rase_minion_rename_reads.py "{input.reads}" \
                        | paste -d '\t' - - \
                        | sort \
                        | uniq \
                        | perl -pe 's/\t/\n/g' \
                        > "{params.reads}"
                else
                    seqtk seq -A "{input.reads}" > "{params.reads}"
                fi
                touch "{output.t}"
            """


rule assign:
    priority: 60
    input:
        ancient("database/.{index}.complete"),
        ancient("matching/.{pref}.read.complete"),
    output:
        t="matching/.{pref}__{index}.assign.complete",
    params:
        reads="matching/{pref}.fa",
        bam="matching/{pref}__{index}.bam",
        index="database/{index}"
    group:
        "group"
    benchmark:
        "benchmarks/{pref}__{index}.classify.log"
    shell:
        """
            prophyle classify -P "{params.index}" -m h1 "{params.reads}" \
                | samtools view -b \
                > "{params.bam}"
            touch "{output.t}"
        """



rule predict:
    priority: 80
    input:
        ancient("matching/.{pref}__{index}.assign.complete"),
        ancient("database/.{index}.complete"),
    output:
        t="prediction/.{pref}__{index}.predict.complete",
    params:
        bam="matching/{pref}__{index}.bam",
        snapshot_dir="prediction/{pref}__{index}/",
        tree="database/{index}/tree.nw",
        metadata="database/{index}.tsv",
        tsv1="prediction/{pref}__{index}.predict_without_flags.tsv",
        tsv2="prediction/{pref}__{index}.predict.tsv",
        reads="matching/{pref}.fa",
    group:
        "group"
    benchmark:
        "benchmarks/{pref}__{index}.predict.log"
    shell:
        """
            mode=read && ../src/rase/rase_has_datetime.py "{params.reads}" || mode=mimic-ont
            mkdir -p "prediction/{wildcards.pref}__{wildcards.index}"
            ../src/rase/rase_predict.py \
                -t $mode \
                -s {config[prediction_sampling]} \
                -p "{params.snapshot_dir}" \
                --pgs-thres-pass {config[pgs_thres_pass]} \
                --ssc-thres-shiconf {config[ssc_thres_Shiconf]} \
                --ssc-thres-sr {config[ssc_thres_SR]} \
                --ssc-thres-rhiconf {config[ssc_thres_Rhiconf]} \
                "{params.tree}" "{params.metadata}" "{params.bam}" \
                > "{params.tsv1}"
            ../src/rase/rase_prediction_add_flags.py  "{params.tsv1}" \
                > "{params.tsv2}"
            rm "{params.tsv1}"
            touch "{output.t}"
        """


rule plot:
    output:
        t="plots/.{pref}__{index}.plot.complete",
    input:
        t1=ancient("plots/.{pref}__{index}.plot_timeline.complete"),
        t2=ancient("plots/.{pref}__{index}.plot_snapshots.complete"),
    group:
        "group"
    shell:
        """
            touch "{output.t}"
        """

rule plot_timeline:
    priority: 90
    output:
        t="plots/.{pref}__{index}.plot_timeline.complete",
    input:
        ancient("prediction/.{pref}__{index}.predict.complete"),
        ancient("database/.{index}.complete"),
    group:
        "group"
    benchmark:
        "benchmarks/{pref}__{index}.plot_timeline.log"
    params:
        tsv="prediction/{pref}__{index}.predict.tsv",
        pdf="plots/{pref}__{index}.timeline.pdf",
        index="database/{index}",
        pref="{pref}"
    shell:
        """
            ../scripts/rase_plot_timeline.R \
                --pgs-thres-pass={config[pgs_thres_pass]} \
                --ssc-thres-shiconf={config[ssc_thres_Shiconf]} \
                --ssc-thres-sr={config[ssc_thres_SR]} \
                --ssc-thres-rhiconf={config[ssc_thres_Rhiconf]} \
                "{params.tsv}" "{params.pdf}"
            touch "{output.t}"
        """


rule plot_snapshots:
    priority: 90
    output:
        t="plots/.{pref}__{index}.plot_snapshots.complete",
    input:
        ancient("prediction/.{pref}__{index}.predict.complete"),
        ancient("database/.{index}.complete"),
    group:
        "group"
    benchmark:
        "benchmarks/{pref}__{index}.plot_snapshots.log"
    params:
        tsv="prediction/{pref}__{index}.predict.tsv",
        index_tsv="database/{index}.tsv",
        pred_dir="prediction/{pref}__{index}",
        plots_pref="plots/{pref}__{index}.snapshots.",
    shell:
        """
            ../src/rase/rase_plot_selected_snapshots.py "{params.index_tsv}" "{params.pred_dir}" 1 5 -1 "{params.plots_pref}"
            touch "{output.t}"
        """


rule decompress:
    priority: 100
    output:
        t="database/.{index}.complete"
    input:
        gz=ancient("database/{index}.tar.gz"),
        tsv=ancient("database/{index}.tsv"),
    benchmark:
        "benchmarks/decompress.{index}.log"
    shell:
        """
            mkdir -p "database/{wildcards.index}"
            prophyle decompress "{input.gz}" database
            touch "{output.t}"
        """
