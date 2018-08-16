.PHONY: all help clean cleanall o2

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p

.SECONDARY:

.SUFFIXES:

all:
	$(SM)

test:
	$(SM) test

replot:
	rm plots/*.pdf
	$(SM)

o2:
	snakemake --cores 9999 -p \
		--cluster-config cluster.o2.json \
		--cluster 'sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem={cluster.memory} --job-name {cluster.name} -o {cluster.output} -e {cluster.error}'

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f plots/*.pdf benchmarks/*.plot.log
	rm -f prediction/.*.quantify.complete benchmarks/*.quantify.log
	rm -f prediction/*.predict.tsv
	rm -f prediction/*.predict.complete
	find prediction -name '*.tsv' | xargs rm -f
	rm -f benchmarks/*.predict.log

cleanall: clean
	rm -f prediction/*.{fq,fq.complete} benchmarks/*.readprep.log
	rm -f prediction/*.bam
	rm -f prediction/.*.complete
	rm -f benchmarks/*.classify.log
	rm -fr database/.*.complete $$(ls -d database/*/ 2>/dev/null || true) benchmarks/decompress.log

