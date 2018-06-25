#! /usr/bin/env python3

import argparse
import collections
import os
import sys
import pysam
import json
import ete3

"""Quantify assignments on the level of leaves.

Assignment: Dictionary with entries name, qlen, h1, c1.
"""

FAKE_ISOLATE_UNASSIGNED="_unassigned_"

def timestamp_from_qname(qname):
    return int(qname.partition("_")[0])


def get_first_timestamp(bam_fn):
    return 0
    f=pysam.AlignmentFile(bam_fn, "rb")
    read = next(f.fetch(until_eof=True))
    f.close()
    return timestamp_from_qname(read.qname)


def format_time(seconds):
    minutes=seconds//60
    hours=minutes//60
    minutes-=60*hours
    return "{}h{}m".format(hours, minutes)


class Stats:
    """Statistics for an experiment.

    Params:
        tree_fn (str): Filename of the Newick tree.

    Attributes:
        nb_reads (int): Number of processed reads.
        nb_nonprop_asgs (int): Number of processed alignments (before propagation).
        nb_asgs (int): Number of processed alignments (after propagation).

        cumul_qlen (float): Cumulative weighted read length.
        cumul_h1 (float): Cumulative weighted hit count.
        cumul_c1 (float): Cumulative weighted coverage.

        stats_count (dict): isolate -> number of processed reads
        stats_qlen (dict): isolate -> weighted qlen
        stats_h1 (dict): isolate -> squared h1
        stats_c1 (dict): isolate -> squared c1
    """

    def __init__(self, tree_fn):
        self.tree=ete3.Tree(tree_fn, format=1)
        self.isolates=sorted([isolate.name for isolate in self.tree])
        self.descending_isolates=self.precompute_descendants(self.tree)

        # stats for assigned reads
        self.nb_assigned_reads=0
        self.nb_nonprop_asgs=0
        self.nb_asgs=0

        # stats for unassigned reads
        self.nb_unassigned_reads=0
       
        # cumulative statistics for assigned reads
        self.cumul_h1_pow1=0.0
        self.cumul_c1_pow1=0.0
        self.cumul_ln_pow1=0.0

        # statistics for individual isolates, "_unassigned_" for unassigned
        self.stats_h1_pow0=collections.defaultdict(lambda:0.0)
        self.stats_h1_pow1=collections.defaultdict(lambda:0.0)
        self.stats_c1_pow1=collections.defaultdict(lambda:0.0)
        self.stats_ln_pow1=collections.defaultdict(lambda:0.0)

    def precompute_descendants(self, tree):
        descending_leaves={}
        for root in list(tree.traverse())+[tree]:
            descending_leaves[root.name]=set([isolate.name for isolate in root])
        return descending_leaves

    def get_number_of_assigned_strains(self, asgs):
        l=0
        for asg in asgs: 
            l+=len(self.descending_isolates[asg["rname"]])
        return l

    def update_from_one_read(self, asgs):
        """Update statistics from assignments of a single read.

        Params:
            asgs (dict): Assignments.
        """

        assert len(asgs)>0, "No assignments provided"

        is_assigned=asgs[0]["assigned"]

        if is_assigned:
            l=self.get_number_of_assigned_strains(asgs)
            
            self.nb_assigned_reads+=1
            self.nb_nonprop_asgs+=len(asgs)
            self.nb_asgs+=l

            for asg in asgs:
                nname=asg["rname"]
                descending_isolates=self.descending_isolates[nname]
                self.update_cumuls(h1=asg["h1"], c1=asg["c1"], ln=asg["ln"], weight=len(descending_isolates)/l)
                self.update_strain_stats(descending_isolates, h1=asg["h1"], c1=asg["c1"], ln=asg["ln"], l=l)

        else:
            assert len(asgs)==1, "A single read shouldn't be reported as unassigned mutliple times (error: {})".format(asgs)
            asg=asgs[0]
            self.nb_unassigned_reads+=1
            self.update_strain_stats([FAKE_ISOLATE_UNASSIGNED], h1=asg["h1"], c1=asg["c1"], ln=asg["ln"], l=l)

    def update_cumuls(self, h1, c1, ln, weight):
        self.cumul_h1_pow1+= h1 * weight
        self.cumul_c1_pow1+= c1 * weight
        self.cumul_ln_pow1+= ln * weight

    def update_strain_stats(self, isolates, h1, c1, ln, l):
        for isolate in isolates:
            self.stats_h1_pow0[isolate]+=1.0 / l
            self.stats_h1_pow1[isolate]+=1.0 * (h1 / l)
            self.stats_c1_pow1[isolate]+=1.0 * (c1 / l)
            self.stats_ln_pow1[isolate]+=1.0 * (ln / l)

    def print(self, file):
        """Print statistics.

        Args:
            file (file): Output file.
        """
        print("taxid", "count", "count_norm", "ln", "ln_norm", "h1", "h1_norm", "c1", "c1_norm", sep="\t", file=file)
        table=[]
        for isolate in self.isolates + [FAKE_ISOLATE_UNASSIGNED]:
            table.append(
                [
                    isolate,
                    self.stats_h1_pow0[isolate],
                    1.0*self.stats_h1_pow0[isolate]/self.nb_assigned_reads if self.nb_assigned_reads!=0 else 0,
                    self.stats_ln_pow1[isolate],
                    self.stats_ln_pow1[isolate]/self.cumul_ln_pow1 if self.cumul_ln_pow1!=0 else 0,
                    self.stats_h1_pow1[isolate],
                    self.stats_h1_pow1[isolate]/self.cumul_h1_pow1 if self.cumul_h1_pow1!=0 else 0,
                    self.stats_c1_pow1[isolate],
                    self.stats_c1_pow1[isolate]/self.cumul_c1_pow1 if self.cumul_c1_pow1!=0 else 0,
                ]
            )

        table.sort(key=lambda x: x[5], reverse=True)

        for x in table:
            print(*x, sep="\t", file=file)


class AssignmentBlockReader:
    """Iterator over blocks of ProPhyle assignments in a BAM file.

    Assumes a non-empty BAM file.

    Params:
        bam_fn (str): BAM file name.

    Attributes:
    """

    def __init__(self, bam_fn):
        self.assignment_reader=AssignmentReader(bam_fn)
        self._buffer=[]
        self._finished=False

    def __iter__(self):
        return self

    def __next__(self):
        """Get next block of assignments of the same read.
        """

        if self._finished:
            raise StopIteration

        while len(self._buffer)<2 or self._buffer[-1]["qname"]==self._buffer[-2]["qname"]:
            try:
                asg=next(self.assignment_reader)
                self._buffer.append(asg)
            except StopIteration:
                self._finished=True
                buffer=self._buffer
                self._buffer=[]
                return buffer

        buffer=self._buffer[:-1]
        self._buffer=self._buffer[-1:]
        return buffer


class AssignmentReader:
    """Iterator over individual ProPhyle assignments in a BAM file.

    Assumes that it is possible to infer read lengths (either
    from the ln tag, base sequence or cigars).

    Params:
        bam_fn (str): BAM file name.

    Attributes:
        samfile (pysam.AlignmentFile): PySAM file object.
    """

    def __init__(self, bam_fn):
        self.samfile = pysam.AlignmentFile(bam_fn, "rb")
        self.alignments_iter = self.samfile.fetch(until_eof=True)

        self.qname=None
        self.last_qname=None
        self.read_ln=None

    def __iter__(self):
        return self

    def __next__(self):
        """Get next assignment.

        Returns:
            assignment (dict): A dict with the following keys: "rname", "qname", "qlen, "assigned", "h1", "c1".
        """

        # 1) acquire a new alignment from SAM
        try:
            alignment=next(self.alignments_iter)
        except StopIteration:
            raise StopIteration

        # 2) re-compute basic read-related variables if necessary
        self.qname=alignment.qname
        if self.qname!=self.last_qname:
            try:
                self.read_ln=alignment.get_tag("ln")
            except KeyError: # a ln tag is not present
                if alignment.seq!="*":
                    self.read_ln=len(alignment.seq)
                else:
                    #self.read_ln=read.infer_read_length()
                    self.read_ln=read.infer_read_length()
        self.last_qname=self.qname

        # 3) infer assignment-related variables
        if not alignment.is_unmapped:
            asg = {
                "rname": alignment.reference_name,
                "qname": alignment.qname,
                "assigned": True,
                "ln": self.read_ln,
                "h1": alignment.get_tag("h1"),
                "c1": alignment.get_tag("c1"),
            }
        else:
            asg = {
                "rname": None,
                "qname": alignment.qname,
                "assigned": False,
                "ln": self.read_ln,
                "h1": None,
                "c1": None,
            }
            print(asg)

        return asg


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
    stats=Stats(args.tree)

    # t=0 point hardcoded as the time of the first read minus 60 seconds
    first_timestamp=get_first_timestamp(args.bam)-60
    print_timestamp=first_timestamp

    # 1) print empty statistics
    f=open("{}/{}.tsv".format(args.pref, print_timestamp), mode="w")

    # 2) iterate through individual reads, and update and print statistics
    assignment_bam_reader=AssignmentBlockReader(args.bam)
    for read_stats in assignment_bam_reader:
        assert len(read_stats)>0
        read_timestamp=timestamp_from_qname(read_stats[0]["qname"])

        if args.pref is not None:

            if print_timestamp + args.delta <= read_timestamp:
                # print current statistics, shift to the new time and open a new file
                stats.print(file=f)
                f.close()
                while print_timestamp + args.delta < read_timestamp:
                    print_timestamp+=args.delta
                f=open("{}/{}.tsv".format(args.pref, print_timestamp), mode="w")

                time=format_time(print_timestamp-first_timestamp)
                print("Time t={}: {} reads and {} non-propagated ({} propagated) assignments processed.".format(time, stats.nb_assigned_reads, stats.nb_nonprop_asgs, stats.nb_asgs), file=sys.stderr)

        stats.update_from_one_read(read_stats)

    stats.print(file=f)
    f.close()

    with open(args.tsv, mode="w") as f:
        stats.print(file=f)

if __name__ == "__main__":
    main()
