#! /usr/bin/env python3

import argparse
import os
import sys
import re
import csv

re_number=re.compile(r'^([0-9]*\.{0,1}[0-9]*)$')
re_greater=re.compile(r'^>={0,1}([0-9]*\.{0,1}[0-9]*)$')
re_lesser=re.compile(r'^<={0,1}([0-9]*\.{0,1}[0-9]*)$')


def pseudo_mic_to_cat(pseudo_mic, breakpoint):
    pseudo_mic=pseudo_mic.strip()
    pseudo_mic=pseudo_mic.replace(" ","")
    pseudo_mic=pseudo_mic.split("/")[0]

    if pseudo_mic=="" or pseudo_mic=="-":
        #return "M"
        return "NA"

    m=re_number.match(pseudo_mic)
    if m:
        number=float(m.group(1))
        if number >= breakpoint:
            return "R"
        else:
            return "S"

    m=re_greater.match(pseudo_mic)
    if m:
        number=float(m.group(1))
        if number >= breakpoint:
            return "R"

    m=re_lesser.match(pseudo_mic)
    if m:
        number=float(m.group(1))
        if number < breakpoint:
            return "S"

    print("MIC value '{}' is ambiguous and resistance could not be resolved".format(pseudo_mic),file=sys.stderr)
    #return "D"
    return "NA"


def assign_cat(
            tsv_in_fn,
            ant,
            breakpoint,
            miccol,
            taxidcol,
        ):

    with open(tsv_in_fn) as tsv_fo:
        #reader = csv.DictReader(tsv_fo, delimiter='\t')

        print("\t".join(["taxid", ant+"_mic", ant+"_cat"]))

        for i,row in enumerate(tsv_fo):
            if i==0:
                continue

            parts=row.split("\t")
            #print(row)
            taxid=parts[taxidcol].strip()
            pseudo_mic=parts[miccol].strip()
            cat=pseudo_mic_to_cat(pseudo_mic, breakpoint)

            if pseudo_mic=="":
                pseudo_mic="NA"

            print("\t".join([taxid, pseudo_mic, cat]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Read MIC's from a TSV file and assign resistance categories (R / S / NA).")

    parser.add_argument('-i', '--in-tsv',
            type=str,
            metavar='str',
            required=True,
            dest='tsv_in_fn',
            help='Input TSV file with taxid / pseudo-mic / cat',
        )

    parser.add_argument('-t', '--taxid-col',
            type=int,
            metavar='int',
            required=True,
            dest='taxidcol',
            help='0-based taxid column ID',
        )

    parser.add_argument('-m', '--mic-col',
            type=int,
            metavar='str',
            required=True,
            dest='miccol',
            help='0-based MIC column ID',
        )

    parser.add_argument('-b', '--breakpoint',
            type=float,
            metavar='float',
            required=True,
            dest='breakpoint',
            help='Breakpoint for antibiotic resistance',
        )

    parser.add_argument('-a', '--antibiotic',
            type=str,
            metavar='str',
            default='ant',
            dest='ant',
            help='Antibiotic acronym (to print) [ant]',
        )

    args = parser.parse_args()

    assign_cat(
            tsv_in_fn=args.tsv_in_fn,
            ant=args.ant,
            breakpoint=args.breakpoint,
            miccol=args.miccol,
            taxidcol=args.taxidcol,
        )
