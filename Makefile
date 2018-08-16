.PHONY: all help clean cleanall o2 export

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p
EXPORTED=export.tar

.SECONDARY:

.SUFFIXES:

all: ## Run everything
	$(SM)

test: ## Run the smallest experiment
	$(SM) test

replot: ## Replot all figures
	rm plots/*.pdf
	$(SM)

export: ## Export all outputs
	rm -f $(EXPORTED)
	tar cvf $(EXPORTED) --files-from /dev/null
	for x in benchmarks/*.log prediction/*.tsv plots/*.pdf; do \
		tar --append --file=$(EXPORTED) $$x; \
	done



o2: ## Submit jobs to Harvard O2
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
	rm -f $(EXPORT)

cleanall: ## Clean all files (including bam files and logs)
cleanall: clean
	rm -f prediction/*.{fq,fq.complete} benchmarks/*.readprep.log
	rm -f prediction/*.bam
	rm -f prediction/.*.complete
	rm -f benchmarks/*.classify.log
	rm -fr database/.*.complete $$(ls -d database/*/ 2>/dev/null || true) benchmarks/decompress.log

