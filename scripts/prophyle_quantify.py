#! /usr/bin/env python3

import argparse
import collections
import os
import sys
import pysam
import json
import ete3

"""Quantify assignments on the level of leaves.

Assignment: Dictionary with entries name, rlen, h1, c1.
"""

class Stats:
    """Statistics for an experiment.

    Params:
        tree_fn (str): Filename of the Newick tree.

    Attributes:
        self.nb_reads (int): Number of processed reads.
        self.nb_nonprop_asgs (int): Number of processed alignments (before propagation).
        self.nb_asgs (int): Number of processed alignments (after propagation).

        self.cumul_rlen (float): Cumulative weighted read length.
        self.cumul_h1 (float): Cumulative weighted hit count.
        self.cumul_c1 (float): Cumulative weighted coverage.

        self.stats_count (dict): nodename -> number of processed reads
        self.stats_rlen (dict): nodename -> weighted rlen
        self.stats_rlensq (dict): nodename -> weighted squared rlen
        self.stats_h1 (dict): nodename -> squared h1
        self.stats_h1sq (dict): nodename -> weighted squared h1
        self.stats_c1 (dict): nodename -> squared c1
        self.stats_c1sq (dict): nodename -> weighted squared c1
    """

    def __init__(self, tree_fn):
        self.tree=ete3.Tree(tree_fn, format=1)

        self.nb_reads=0
        self.nb_nonprop_asgs=0
        self.nb_asgs=0

        self.cumul_rlen=0.0
        self.cumul_h1=0.0
        self.cumul_c1=0.0

        self.stats_count=collections.defaultdict(lambda:0.0)
        self.stats_rlen=collections.defaultdict(lambda:0.0)
        self.stats_rlensq=collections.defaultdict(lambda:0.0)
        self.stats_h1=collections.defaultdict(lambda:0.0)
        self.stats_h1sq=collections.defaultdict(lambda:0.0)
        self.stats_c1=collections.defaultdict(lambda:0.0)
        self.stats_c1sq=collections.defaultdict(lambda:0.0)

        self.nodename_to_leave_nodenames=collections.defaultdict(lambda:[])
        for root in list(self.tree.traverse())+[self.tree]:
            self.nodename_to_leave_nodenames[root.name]=set([x.name for x in root])

        #print(self.nodename_to_leave_nodenames, file=sys.stderr)


    def propagate_asgs_leaves(self, asgs):
        """Propagate assignments of a single from nodes to leaves.

        Params:
            asgs (list of dicts): Assignments.

        Return:
            asgs_leaves (list of dicts): Assignments propagated to the leaves.
        """

        asgs_leaves=[]
        for asg in asgs:
            for leafname in self.nodename_to_leave_nodenames[asg['name']]:
                    asg2=asg.copy()
                    asg2["name"]=leafname
                    asgs_leaves.append(asg2)

        return asgs_leaves


    def update_oneread(self, asgs):
        """Update statistics from assignments of a single read.

        Params:
            asgs (dict): Assignments.
        """

        self.nb_nonprop_asgs+=len(asgs)
        asgs_leaves=self.propagate_asgs_leaves(asgs)
        l=len(asgs_leaves)

        self.nb_reads+=1
        self.nb_asgs+=l

        # in theory, h1 and c1 can be different for different assignments
        self.cumul_rlen+=1.0*sum([x["rlen"] for x in asgs_leaves]) / l
        self.cumul_h1+=1.0*sum([x["h1"] for x in asgs_leaves]) / l
        self.cumul_c1+=1.0*sum([x["c1"] for x in asgs_leaves]) / l

        for asg in asgs_leaves:
            n=asg["name"]
            self.stats_rlen[n]+=1.0 * asg["rlen"] / l
            self.stats_rlensq[n]+=1.0 * (asg["rlen"]**2) / l
            self.stats_count[n]+=1.0 / l
            self.stats_h1[n]+=1.0*asg["h1"] / l
            self.stats_h1sq[n]+=1.0*(asg["h1"]**2) / l
            self.stats_c1[n]+=1.0*asg["c1"] / l
            self.stats_c1sq[n]+=1.0*(asg["c1"]**2) / l


    def print(self, file):
        """Print statistics.

        Args:
            file (file): Output file.
        """
        print("taxid", "count", "count_norm", "len", "len_norm", "lensq", "h1", "h1_norm", "h1sq", "c1", "c1_norm", "c1sq", sep="\t", file=file)
        table=[]
        for node in self.tree:
            n=node.name
            table.append(
                [
                    n,
                    self.stats_count[n],
                    1.0*self.stats_count[n]/self.nb_reads if self.nb_reads!=0 else 0,
                    self.stats_rlen[n],
                    self.stats_rlen[n]/self.cumul_rlen if self.cumul_rlen!=0 else 0,
                    self.stats_rlensq[n],
                    self.stats_h1[n],
                    self.stats_h1[n]/self.cumul_h1 if self.cumul_h1!=0 else 0,
                    self.stats_h1sq[n],
                    self.stats_c1[n],
                    self.stats_c1[n]/self.cumul_c1 if self.cumul_c1!=0 else 0,
                    self.stats_c1sq[n],
                ]
            )

        table.sort(key=lambda x: x[7], reverse=True)

        for x in table:
            print(*x, sep="\t", file=file)


class BamReader:
    """BAM Reader.

    Params:
        bam_fn (str): BAM file name.

    Attributes:
        samfile (pysam.AlignmentFile): PySAM file object.
        asgs (list of dicts): Assignments (unfinished).
    """

    def __init__(self, bam_fn):
        self.samfile = pysam.AlignmentFile(bam_fn, "rb")
        self.asgs=collections.defaultdict(lambda:[])


    def process_read(self):
        """Generator function returning assignments of a single reads.

        Returns:
            assignments (list of dicts): Assignments (all levels, non-propagated).
        """

        self.name=None
        last_name=None
        read_ln=0
        asgs=[]

        for read in self.samfile.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            self.name=read.qname

            if self.name!=last_name:
                read_ln=read.infer_read_length()

            asg = {
                "name": read.reference_name,
                "rlen": read_ln,
                "h1": read.get_tag("h1"),
                "c1": read.get_tag("c1"),
            }

            if last_name is None or self.name==last_name:
                asgs.append(asg)
            else:
                #print(asgs)
                yield asgs
                asgs=[asg]
            last_name=self.name
        return asgs

    def __del__(self):
        self.samfile.close()


def timestamp_from_qname(qname):
    return int(qname.partition("_")[0])


def get_first_timestamp(bam_fn):
    return 0
    bamreader=BamReader(bam_fn)
    for read_stats in bamreader.process_read():
        return timestamp_from_qname(bamreader.name)


def format_time(seconds):
    minutes=seconds//60
    hours=minutes//60
    minutes-=60*hours
    return "{}h{}m".format(hours, minutes)


def main():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument('tree',
            type=str,
            metavar='<tree.nw>',
        )

    parser.add_argument('bam',
            type=argparse.FileType('r'),
            metavar='<assignments.bam>',
        )

    parser.add_argument('tsv',
            type=str,
            metavar='stats.tsv',
        )

    parser.add_argument('-p',
            type=str,
            dest='pref',
            metavar='STR',
            default=None,
            help="Output dir for samplings"
        )

    parser.add_argument('-i',
            type=int,
            dest='delta',
            metavar='INT',
            help='sampling interval (in seconds) [300]',
            default=300,
        )


    args = parser.parse_args()
    bamreader=BamReader(args.bam)
    stats=Stats(args.tree)

    # t=0 point hardcoded as the time of the first read minus 60 seconds
    first_timestamp=get_first_timestamp(args.bam)-60
    print_timestamp=first_timestamp

    # 1) print empty statistics
    f=open("{}/{}.tsv".format(args.pref, print_timestamp), mode="w")

    # 2) iterate through individual reads, and update and print statistics
    for read_stats in bamreader.process_read():
        read_timestamp=timestamp_from_qname(bamreader.name)

        if args.pref is not None:

            if print_timestamp + args.delta <= read_timestamp:
                # print current statistics, shift to the new time and open a new file
                stats.print(file=f)
                f.close()
                while print_timestamp + args.delta < read_timestamp:
                    print_timestamp+=args.delta
                f=open("{}/{}.tsv".format(args.pref, print_timestamp), mode="w")

                time=format_time(print_timestamp-first_timestamp)
                print("Time t={}: {} reads and {} non-propagated ({} propagated) assignments processed.".format(time, stats.nb_reads, stats.nb_nonprop_asgs, stats.nb_asgs), file=sys.stderr)

        stats.update_oneread(read_stats)

    stats.print(file=f)
    f.close()

    with open(args.tsv, mode="w") as f:
        stats.print(file=f)

if __name__ == "__main__":
    main()
