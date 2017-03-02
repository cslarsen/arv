# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

import arv
import unittest

class ArvModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("fake_genome.txt")

    def test_module_len(self):
        self.assertEqual(len(self.genome), 13)
        self.assertEqual(len(arv.Genome()), 0)

    def test_module_ychromo(self):
        self.assertTrue(self.genome.ychromo)
        self.assertFalse(arv.Genome().ychromo)

    def test_module_genotypes(self):
        self.assertEqual(self.genome[4477212], "AT")
        self.assertEqual(self.genome["rs4477212"], "AT")
        self.assertEqual(self.genome["rs4672279"], "GT")
        self.assertEqual(self.genome["rs742927"], "GG")

        with self.assertRaises(TypeError):
            self.genome[1.0]

        with self.assertRaises(KeyError):
            self.genome[123]

        with self.assertRaises(KeyError):
            self.genome["rs123"]

    def test_module_load_factor(self):
        self.assertIsInstance(self.genome.load_factor(), float)
        self.assertGreater(self.genome.load_factor(), 0.0)
        self.assertLess(self.genome.load_factor(), 1.0)
