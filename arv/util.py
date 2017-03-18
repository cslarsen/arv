# Copyright 2014, 2016, 2017 Christian Stigen Larsen
# Distributed under the GPL v3 or later. See COPYING.

def make_report(genome, functions):
    """Runs each function with genome as argument, returning a dict of
    results."""
    report = {}

    for func in functions:
        if func.__doc__ is not None:
            title = func.__doc__[:func.__doc__.index(".")]
        else:
            title = func.__name__.replace("_", " ").capitalize()

        try:
            report[title] = func(genome)
        except ValueError as e:
            report[title] = "Error: %s" % e
        except AssertionError as e:
            report[title] = "Error: %s" % e
        except NotImplementedError:
            continue

    return report
