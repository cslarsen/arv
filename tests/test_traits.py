"""
Inferring tests for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import arv
import arv.traits
import unittest

class ArvTraitsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("tests/fake_genome.txt", ethnicity="european")

    def test_alcohol_flush_reaction(self):
        self.assertEqual(self.genome["rs671"], "GG")
        self.assertEqual(arv.traits.alcohol_flush_reaction(self.genome),
                "Little to no reaction (two copies of the ALDH2 gene)")
