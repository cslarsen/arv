"""
Various benchmarks for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import argparse
import arv
import contextlib
import os
import random
import sys
import time
import unittest

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

def log(msg, stream=sys.stdout):
    stream.write(msg)
    stream.flush()

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
    stream = local_args.get("stream", sys.stdout)
    prefix = local_args.get("prefix", "")
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
            if round(elapsed, 4) < round(best, 4):
                log("\n%s%6.4fs " % (prefix, elapsed), stream=stream)
            best = elapsed
        else:
            log(".", stream=stream)
    return best

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

class BenchmarkTests(unittest.TestCase):
    @unittest.skipUnless(os.path.isfile(os.getenv("ARV_BENCHMARK", ".")),
        "Specify ARV_BENCHMARK=<genome-filename> to benchmark")
    def test_parser_speed(self):
        times = 40
        code = benchmarks["parsing"]
        filename = os.getenv("ARV_BENCHMARK")
        seconds = benchmark(times, code, filename=filename, stream=sys.stderr,
                prefix="  ")
        genome = arv.load(filename)
        sys.stderr.flush()
        sys.stderr.write(" %d SNPs in ~%dms or %.1g SNPs/second ... " % (
                len(genome), int(round(seconds, 3)*1000), len(genome)/seconds))
        sys.stderr.flush()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--times", "-t", default=20, type=int)
    p.add_argument("--filename", "-f", required=True, type=str)
    args = p.parse_args()
    all_benchmarks(args.filename, args.times)
