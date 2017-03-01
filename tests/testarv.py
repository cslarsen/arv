# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

import arv
import unittest

class ArvTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("fake_genome.txt")

    def test_len(self):
        self.assertEqual(len(self.genome), 13)
        self.assertEqual(len(arv.PyGenome()), 0)

    def test_ychromo(self):
        self.assertTrue(self.genome.ychromo)
        self.assertFalse(arv.PyGenome().ychromo)

    def test_genotypes(self):
        self.assertEqual(self.genome[4477212], "AT")
        self.assertEqual(self.genome["rs4477212"], "AT")
        self.assertEqual(self.genome["rs4672279"], "GT")
        self.assertEqual(self.genome["rs742927"], "GG")

