#! /usr/bin/env python3

import argparse
import datetime
import json
import os
import re
import sys

#rre=re.compile(r'^>([^\s]+).*start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})')
#rre=re.compile(r'^>([0-9a-f\-]+) runid=([0-9a-f]) read=(\d+) ch=(\d+) start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})')

#rre=re.compile(r'''^>([0-9a-f\-]+) runid=([0-9a-f]+) read=(\d+) ch=(\d+) start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})''')#, flags=re.X)
#rre=re.compile(r'''^>([^\s]+)(.*) \s
#    start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})''', flags=re.X)

re_minion=re.compile(r'''^[>@]
        ([0-9a-f\-]+) \s
        runid=([0-9a-f]+) \s
        read=(\d+) \s
        ch=(\d+) \s
        start_time=(\d{4}-\d{2}-\d{2})T(\d{2}:\d{2}:\d{2})
    ''', flags=re.X)



def parse_header_line(hl):
    """Take a FA/FQ header line and parse all Minion fields.

    Args:
        hl (str): Header line.

    Returns:
        (uuid, runid, read, ch, date, time)
    """
    m=re_minion.match(hl)
    gg=m.groups()
    assert gg is not None, "Header line doesn't match the regular expression ('{}').".format(hl)
    return gg


def timestamp_from_datetime(date, time):
    """Get a timestamp from a date and time.

    Args:
        date (str): Date ('YY-MM-DD').
        time (str): Time ('hh:mm:ss').

    Returns:
        timestamp (int): Unix timestamp.
    """
    date_time=date+" "+time
    b=datetime.datetime.strptime(date_time, "%Y-%m-%d %H:%M:%S" )
    ts=int(b.timestamp())
    return ts


def transf(fq):
    d={}

    with open(fq) as f:
        for i,x in enumerate(f):
            if i % 4 == 0:
                uuid, runid, read, ch, date, time = parse_header_line(x)
                ts=timestamp_from_datetime(date, time)
                print("@"+"{}_ex{}_r{}_ch{}".format(ts, runid[:9], uuid.partition("-")[0], ch))
            else:
                print(x,end="")



def main():
    parser = argparse.ArgumentParser(description="Reformat Minion read names in a FQ file.")

    parser.add_argument('fq',
            type=str,
            metavar='reads.fq',
        )

    args = parser.parse_args()

    transf(args.fq)

if __name__ == "__main__":
    main()
