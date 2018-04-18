#! /usr/bin/env python3

import argparse
import collections
import csv
import itertools
import os
import re
import sys

re_cat=re.compile("(.*)_cat")
cats=['S', 'R', 'NA']

def main():

    stats=collections.defaultdict(lambda: collections.Counter())
    ants=collections.OrderedDict()
    pgs=set()
    pg_count=collections.Counter()

    with open("res_cat.tsv") as f:
        tsv_reader = csv.DictReader(f, delimiter='\t')
        for r in tsv_reader:
            pg=r["phylogroup"]
            pg_count[pg]+=1
            pgs.add(pg)

            for k in r:
                m=re_cat.match(k)
                if m:
                    ant=m.group(1)
                    cat=r[k]
                    ants[ant]=''

                    stats[pg][(ant, cat)]+=1

            #print(pg,r)

    parts=['pg', 'count'] + ['{}_{}'.format(ant, cat) for ant in ants for cat in cats]
    print(*parts, sep="\t")

    for pg in sorted(pgs, key=lambda x: int(x)):
        parts=[pg, pg_count[pg]]
        for ant in ants:
            for cat in cats:
                parts.append(stats[pg][(ant,cat)])

        print(*parts, sep="\t")

if __name__ == "__main__":
    main()
