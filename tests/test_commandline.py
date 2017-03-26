"""
Tests for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import arv
import os
import subprocess
import sys
import unittest

class ArvModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome_path = os.path.join(os.path.dirname(__file__),
                "fake_genome.txt")

    def _execute(self, *args):
        output = subprocess.check_output([sys.executable, "-m", "arv"] +
                list(args), universal_newlines=True)
        return output.replace("\r\n", "\n").split("\n")

    def test_help(self):
        self.assertTrue("\n".join(self._execute("--help")).
                startswith("usage: arv [-h]"))

    def test_example(self):
        self.assertEqual(self._execute("--example", self.genome_path),
            ["fake_genome.txt ... 25 SNPs, male",
             "fake_genome.txt ... A man with blue eyes and light skin",
             ""])
