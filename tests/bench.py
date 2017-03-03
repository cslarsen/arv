#! /usr/bin/env python
# -*- encoding: utf-8 -*-

"""
Benchmarking tool for arv

Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later. See COPYING for details.
"""

import argparse
import os
import arv
import sys
import time

def parse_args():
    p = argparse.ArgumentParser(description="arv parse benchmarking")
    p.add_argument("--count", "-c", type=int, metavar="N", default=80,
        help="Exit when time hasn't improved after N tries")
    p.add_argument("file", type=str, help="Genome file to parse")

    opts = p.parse_args(sys.argv[1:])

    if not os.path.isfile(opts.file):
        print("No such file: %s" % opts.file)
        sys.exit(1)

    if opts.count < 0:
        print("Invalid --count value: %d" % opts.count)
        sys.exit(1)

    return opts

def benchmark(filename, maxsteps=20):
    best = 9999
    steps = 0
    while True:
        start = time.clock()
        genome = arv.load(filename)
        stop = time.clock()

        this = stop-start
        if this < best:
            best = this
            sys.stdout.write("\n%.4fs " % (stop-start))
            sys.stdout.flush()
            steps = 0
        else:
            steps += 1
            sys.stdout.write(".")
            sys.stdout.flush()
            if steps >= maxsteps:
                break
    print("")
    return best

if __name__ == "__main__":
    opts = parse_args()
    sys.stdout.write("Benchmarking arv.load(%r) " % opts.file)
    best = benchmark(opts.file, opts.count)
    print("Best time: %fs" % best)
