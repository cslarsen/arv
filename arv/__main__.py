"""
Command line interface to arv.

This can be invoked with ``python -m arv``.

Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later. See COPYING.
"""

import __init__ as arv
import argparse
import os
import sys

class ArvError(RuntimeError):
    pass

def log(msg="\n"):
    sys.stdout.write(msg)
    sys.stdout.flush()

def _parse_args():
    p = argparse.ArgumentParser(prog="arv",
            description="Arv - a fast 23andMe parser",
            epilog=arv.__copyright__)

    p.add_argument("--repl", default=False, action="store_true",
            help="Open a Python REPL loaded with the given genomes")

    p.add_argument("--example", default=False, action="store_true",
            help="Shows an example report for the genome(s)")

    p.add_argument("files", nargs="+",
            help="23andMe raw genome file name(s)")

    opts = p.parse_args()
    return opts

def summary(genome):
    """Returns a textual summary of the genome."""
    return "{count} SNPs, {gender}".format(
            count=len(genome), gender="male" if
            genome.y_chromosome else "female",)

def example(genome):
    """Returns an example report for the genome."""
    gender = "man" if genome.y_chromosome else "woman"
    complexion = "light" if genome["rs1426654"] == "AA" else "dark"

    color = arv.unphased_match(genome["rs12913832"], {
        "AA": "brown",
        "AG": "brown or green",
        "GG": "blue"})

    report = "A {gender} with {color} eyes and {complexion} skin"
    return report.format(**locals())

def _main():
    opts = _parse_args()

    genomes = []
    for filename in opts.files:
        log("%s ... " % os.path.basename(filename))
        genome = arv.load(filename)
        log("%s\n" % summary(genome))
        genomes.append(genome)

    if opts.example:
        log()
        for filename, genome in zip(opts.files, genomes):
            log("%s ... %s\n" % (os.path.basename(filename), example(genome)))

    if opts.repl:
        env = dict(globals())

        if len(genomes) == 1:
            env.update({"genome": genomes[0]})
            message = "Type `genome` to see the parsed 23andMe raw genome file"
        else:
            env.update({"genomes": genomes})
            message = "Type `genomes` to see the parsed 23andMe raw genome files"

        import code
        code.interact(message, local=env)

if __name__ == "__main__":
    try:
        _main()
        sys.exit(0)
    except ArvError as e:
        log("Error: %s\n" % e)
        sys.exit(1)
