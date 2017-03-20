"""
Various benchmarks for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import argparse
import arv
import contextlib
import random
import sys
import time

benchmarks = {
    "parsing": "arv.load(filename)",

    "random access":
r"""
for n in xrange(5000):
    try:
        pos = random.randint(genome.first, genome.last)
        snp = genome[pos]
    except KeyError:
        # RSIDs are not contiguous
        pass
""",

    "iterate items in genome":
r"""
assert(False) # this is too slow at the moment
num = 0
for snp in genome:
    num += 1
assert(num == len(genome))
""",

    "iterate rsids":
r"""
num = 0
for snp in genome.rsids:
    num += 1
assert(num == len(genome))
""",

    "iterate snps":
r"""
num = 0
for snp in genome.snps:
    num += 1
assert(num == len(genome))
""",
}

def log(msg):
    sys.stdout.write(msg)
    sys.stdout.flush()

if sys.version_info[:2] >= (3, 3):
    mark_time = time.perf_counter
else:
    mark_time = time.clock

@contextlib.contextmanager
def timed_block():
    start = mark_time()
    elapsed = None
    yield lambda: elapsed
    elapsed = mark_time() - start

def benchmark(times, code, **local_args):
    best = 1e9
    for no in range(times):
        localvars = {
            "arv": arv,
            "sys": sys
        }
        if sys.version_info[0] >= 3:
            localvars["xrange"] = range

        localvars.update(local_args)

        with timed_block() as elapsed:
            exec(code, localvars)

        elapsed = elapsed()
        if elapsed < best:
            log("\n%6.4fs " % elapsed)
            best = elapsed
        else:
            log(".")
    return elapsed

def all_benchmarks(filename, times):
    log("Benchmarking arv %s at %s\n" % (arv.__version__, arv.__file__))
    log("Measuring time with %s\n\n" % mark_time)

    results = {}
    genome = arv.load(filename)

    for name, code in sorted(benchmarks.items()):
        log("Benchmarking %s x %d ... " % (repr(name), times))
        try:
            results[name] = benchmark(times, code, filename=filename,
                    genome=genome, random=random)
        except Exception as e:
            log(str(e))
        finally:
            if name == "parsing":
                genome = arv.load(filename)
                log(" %.2g SNPs / second" % (len(genome)/results[name]))
            log("\n")

    return results

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--times", "-t", default=20, type=int)
    p.add_argument("--filename", "-f", required=True, type=str)
    args = p.parse_args()
    all_benchmarks(args.filename, args.times)
