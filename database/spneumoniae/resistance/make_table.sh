#! /usr/bin/env bash

set -e
set -o pipefail
set -u

set -x

cut -f3 isolates.2013.tsv > _taxids1.txt
cut -f2 isolates.2015.tsv > _taxids2.txt
diff _taxids1.txt _taxids2.txt

./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 0.06 -m 12 -a pnc > _res_cat.pnc.tsv # Benzylpenicillin
./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 0.25 -m 13 -a cft > _res_cat.cft.tsv # Ceftriaxone
./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 1.00 -m 14 -a trm > _res_cat.trm.tsv # Trimethoprim
./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 0.25 -m 15 -a ert > _res_cat.ert.tsv # Erythromycin
./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 1.00 -m 16 -a ttr > _res_cat.ttr.tsv # Tetracycline
./assign_re_cat.py -i isolates.2013.tsv -t 2 -b 8.00 -m 17 -a chl > _res_cat.chl.tsv # Chloramphenicol

for x in _res*.tsv; do
    ./infer_missing_res_cat.py -i $x -t SPARC.core_genes.tree > _$x
done


paste \
    <(cat  _taxids1.txt | sed 's/Taxon ID/taxid/g') \
    <(cut -f 13 isolates.2015.tsv | sed 's/Sequence Cluster/phylogroup/g') \
    <(cut -f 9 isolates.2015.tsv | sed 's/Serotype From Reads/serotype/g') \
    <(cut -f 11 isolates.2015.tsv | sed 's/Sequence Type From Reads/ST/g') \
    <(cut -f 2-3 __res_cat.pnc.tsv) \
    <(cut -f 2-3 __res_cat.cft.tsv) \
    <(cut -f 2-3 __res_cat.trm.tsv) \
    <(cut -f 2-3 __res_cat.ert.tsv) \
    <(cut -f 2-3 __res_cat.ttr.tsv) \
    <(cut -f 2-3 __res_cat.chl.tsv) \
> res_cat.tsv


