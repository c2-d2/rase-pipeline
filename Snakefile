import glob
import os
import sys

snakemake.shell.prefix("set -eo pipefail;")

index="sparc1_k18"

experiments=[os.path.basename(x[:-3]) for x in glob.glob("reads/*.fq")]
print("Experiments:", experiments)


rule all:
    input:
        [
            [
                "reads/preprocessed/{}.fq.complete".format(e),
                "prediction/{}.bam.complete".format(e),
                "prediction/{}.prediction.complete".format(e)
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
        t="reads/preprocessed/{pref}.fq.complete"
    params:
        fq="reads/preprocessed/{pref}.fq"
    shell:
        """
            ./scripts/minion_rename_reads.py {input.fq} | paste -d '\t' - - - - | sort | perl -pe 's/\t/\n/g' > "{params.fq}"
            touch "{output.t}"
        """


rule classify:
    input:
        t="reads/preprocessed/{pref}.fq.complete"
    output:
        t="prediction/{pref}.bam.complete",
    params:
        fq="reads/preprocessed/{pref}.fq",
        bam="prediction/{pref}.bam",
        index="database/{}".format(index)
    shell:
        """
            prophyle classify -P "{params.index}" -m h1 "{params.fq}" \
                | samtools view -b > "{params.bam}"

            touch "{output.t}"
        """


rule quantify_complete:
    input:
        t="prediction/{pref}.bam.complete"
    output:
        t="prediction/{pref}.quantify.complete"
    params:
        pref="{pref}",
        bam="prediction/{pref}.bam",
        index=index
    #benchmark:
        #"benchmarks/quantify__{pref}.{rest}.log"
    shell:
        """
            mkdir -p "prediction/{params.pref}"

            scripts/prophyle_quantify.py -p 'prediction/{params.pref}/' -i 60 "database/{params.index}/tree.nw" "{params.bam}" /dev/null

            touch "{output.t}"
        """


rule predict:
    input:
        t="prediction/{pref}.bam.complete",
    output:
        t="prediction/{pref}.predict.complete",
    benchmark:
        "benchmarks/snapshots__{pref}.{rest}.log"
    shell:
        """
        """


#rule database:
#    input:
#        t="database/.complete"


rule decompress:
    output:
        t="database/{}.complete".format(index)
    input:
        gz="database/{}.tar.gz".format(index),
        tsv="database/{}.tsv".format(index)
    params:
        index=index
    shell:
        """
            prophyle decompress {input.gz} database/{params.index}
            touch "{output.t}"
        """

