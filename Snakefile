import glob
import os
import sys


snakemake.shell.prefix("set -eo pipefail;")


index_tsv=list(map(os.path.basename, glob.glob("database/*.tsv")))
index_tar=list(map(os.path.basename, glob.glob("database/*.tar.gz")))
assert len(index_tar)==1, "The database directory should contain exactly 1 TAR.GZ file (with a ProPhyle compressed index)"
assert len(index_tsv)==1, "The database directory should contain exactly 1 TSV file (with the metadata)"
index_tar=index_tar[0]
index_tsv=index_tsv[0]
index=index_tar.replace(".tar.gz", "")

experiments=[os.path.basename(x[:-3]) for x in glob.glob("reads/*.fq")]
print("Experiments:", experiments)
if len(experiments)==0:
    print("!!!! ", file=sys.stderr)
    print("!!!! WARNING ", file=sys.stderr)
    print("!!!! ", file=sys.stderr)
    print("!!!! No FASTQ files provided", file=sys.stderr)
    print("!!!! ", file=sys.stderr)


rule all:
    input:
        [
            [
                "prediction/{}.fq.complete".format(e),
                "prediction/{}.bam.complete".format(e),
                "prediction/{}.quantify.complete".format(e),
                "prediction/{}.predict.complete".format(e),
                "plots/{}.timeline.pdf".format(e),
            ]
           for e in experiments
        ],
        ##"database/{}.complete".format(index)
        #"prediction/{pref}.bam.complete"
    run:
        print("RASE is starting")


rule preprocess_reads:
    input:
        fq="reads/{pref}.fq"
    output:
        t="prediction/{pref}.fq.complete"
    params:
        fq="prediction/{pref}.fq"
    benchmark:
        "benchmarks/{pref}.readprep.log"
    shell:
        """
            ./scripts/minion_rename_reads.py {input.fq} | paste -d '\t' - - - - | sort | perl -pe 's/\t/\n/g' > "{params.fq}"
            touch "{output.t}"
        """


rule classify:
    input:
        "database/{}.complete".format(index),
        "prediction/{pref}.fq.complete"
    output:
        t="prediction/{pref}.bam.complete",
    params:
        fq="prediction/{pref}.fq",
        bam="prediction/{pref}.bam",
        index="database/{}".format(index)
    benchmark:
        "benchmarks/{pref}.classify.log"
    shell:
        """
            prophyle classify -P "{params.index}" -m h1 "{params.fq}" \
                | samtools view -b > "{params.bam}"

            touch "{output.t}"
        """


rule quantify_complete:
    input:
        "prediction/{pref}.bam.complete",
        "database/{}.complete".format(index),
    output:
        t="prediction/{pref}.quantify.complete"
    params:
        pref="{pref}",
        bam="prediction/{pref}.bam",
        index=index
    benchmark:
        "benchmarks/{pref}.quantify.log"
    shell:
        """
            mkdir -p "prediction/{params.pref}"
            scripts/prophyle_quantify.py -p 'prediction/{params.pref}/' -i 60 "database/{params.index}/tree.nw" "{params.bam}" /dev/null
            touch "{output.t}"
        """


rule predict:
    input:
        "prediction/{pref}.quantify.complete",
        "database/{}.complete".format(index),
    output:
        t="prediction/{pref}.predict.complete",
    params:
        tsv="prediction/{pref}.predict.tsv",
        index=index
    benchmark:
        "benchmarks/{pref}.predict.log"
    shell:
        """
            scripts/rase_predict.py database/{params.index}.tsv prediction/{wildcards.pref}/*.tsv > "{params.tsv}"
            touch "{output.t}"
        """


rule plot_timeline:
    input:
        "prediction/{pref}.predict.complete",
        "database/{}.complete".format(index)
    output:
        pdf="plots/{pref}.timeline.pdf",
    benchmark:
        "benchmarks/{pref}.plot.log"
    params:
        tsv="prediction/{pref}.predict.tsv",
        index="database/{}".format(index),
        pref="{pref}"
    shell:
        """
            ./scripts/plot_selected_snapshots.py {params.index}.tsv prediction/{params.pref} 0 4 -1 plots/{params.pref}.snapshots.
            ./scripts/plot_timeline.R {params.tsv} {output.pdf}
        """


"""
Todo: add possible checking of index consistency
"""
rule decompress:
    output:
        t="database/{}.complete".format(index)
    input:
        gz="database/{}".format(index_tar),
        tsv="database/{}".format(index_tsv)
    benchmark:
        "benchmarks/decompress.log"
    shell:
        """
            prophyle decompress {input.gz} database
            touch "{output.t}"
        """

