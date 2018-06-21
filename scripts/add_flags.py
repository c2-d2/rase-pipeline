#! /usr/bin/env python3

import argparse
import collections
import csv
import os
import re
import sys

# regular expressions for categories to be flagged
res_flagging=[
            re.compile(r"^taxid$"),
            re.compile(r"^serotype$"),
            re.compile(r"^ST$"),
            re.compile(r"^PG\d+$"),
            re.compile(r"_cat$"),
        ]


def add_flags(tsv_fn):
    #load the last record
    with open(tsv_fn) as f:
        for last_rec in csv.DictReader(f, delimiter="\t"):
            pass

    #add flags
    with open(tsv_fn) as f:
        prev_rec=collections.defaultdict(lambda: "")
        for i, rec in enumerate(csv.DictReader(f, delimiter="\t")):
            if i==0:
                print(*rec.keys(), "flags", sep="\t")
                flag_cols=[]
                for k in rec.keys():
                    for re_flagging in res_flagging:
                        if re_flagging.match(k):
                            flag_cols.append(k)
                            break
            flags=[]
            for k in flag_cols:
                if rec[k]!=prev_rec[k]:
                    if rec[k]==last_rec[k]:
                        # detected
                        flags.append(k+'_D')
                    elif rec[k]==prev_rec[k]:
                        # lost successfuly detected
                        flags.append(k+'_L')
            print(*rec.values(), flags, sep="\t")
            prev_rec=rec


def main():
    parser = argparse.ArgumentParser(description="Add flags to predictions (R=beginning of a run, S=stabilized).")

    parser.add_argument('file',
           type=str,
           metavar='predictions.tsv',
           help='',
           )

    args = parser.parse_args()
    add_flags(args.file)

if __name__ == "__main__":
    main()
