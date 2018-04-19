#! /usr/bin/env python3

import argparse
import collections
import datetime
import glob
import os
import re
import sys

re_fn=re.compile(r'(.*?)\..*\.(\d+)\.pdf')

def get_exp_inf(fn):
    """Get experiment identifier.
    """
    bn=fn.split("/")[-1]
    m=re_fn.match(bn)
    exp,ts=m.groups()
    return exp, int(ts)


def seconds_to_text(cs):
    v=collections.OrderedDict()
    v['s']=cs%60
    cs//=60
    v['m']=cs%60
    cs//=60
    v['h']=cs%24
    cs//=24
    v['d']=cs
    parts=[]
    for k in v:
        if v[k]!=0:
            parts.append("{}{}".format(v[k],k))
    return "_".join(parts[::-1])



def plot_snapshots(directory, indexes, outprefix):
    exps=collections.defaultdict(list)

    fns=glob.glob("../*.pdf")
    for fn in fns:
        #print(fn)
        exp,ts=get_exp_inf(fn)
        exps[exp].append([fn, ts])

    for k in exps:
        exps[k].sort(key=lambda x:x[0])

        time0=exps[k][0][1]-first_ts

        for exp, ts in exps[k]:
            abs_ts=ts-time0
            t=seconds_to_text(abs_ts)
            #print(exp, k, t)
            new_fn="{}.{}.pdf".format(k, t)
            print("{} -> {}".format(exp, new_fn))
            shutil.copy(exp,new_fn)


def main():
    parser = argparse.ArgumentParser(description="Plot selected snapshots and use human readable names (e.g., 5m, 1h, etc.)")

    parser.add_argument('-f',
            type=int,
            default=60,
            metavar='int',
            dest='first',
            help='First timestamp [60]',
        )

    parser.add_argument('dir',
            type=str,
            metavar='dir',
            help='Snapshot directory',
        )

    parser.add_argument('indexes',
            type=str,
            metavar='index',
            nargs='+',
            help='0-based index of a snapshot (e.g., 0, -2, 6)',
        )

    parser.add_argument('outpref',
            type=str,
            metavar='outpref',
            help='Output prefix',
        )


    args = parser.parse_args()

    plot_snapshots(args.dir, args.indexes, args.outpref)


if __name__ == "__main__":
    main()

