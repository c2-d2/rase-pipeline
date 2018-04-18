import glob
import os
import sys

snakemake.shell.prefix("set -eo pipefail;")

experiments=[os.path.basename(x[:-3]) for x in glob.glob("reads/*.fq")]
print("Experiments:", experiments)

rule all:
    input:
        ["reads/preprocessed/{}.fq.complete".format(e) for e in experiments]
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

