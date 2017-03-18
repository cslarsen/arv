# arv
# Copyright 2017 Christian Stigen Larsen
# Distributed under the GNU GPL v3 or later; see COPYING.

import arv
import unittest

class ArvModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("tests/fake_genome.txt")

    def test_module_len(self):
        self.assertEqual(len(self.genome), 17)
        self.assertEqual(len(arv.Genome()), 0)

    def test_module_ychromosome(self):
        self.assertTrue(self.genome.y_chromosome)
        self.assertFalse(arv.Genome().y_chromosome)

    def test_module_genotypes(self):
        # Internal IDs not supported yet
        # self.assertEqual(self.genome["i3001754"], "A")
        # self.assertEqual(self.genome["i3001755"], "--")
        # self.assertEqual(self.genome["i3001759"], "--")
        # self.assertEqual(self.genome["i3001761"], "--")
        # self.assertEqual(self.genome["i3001773"], "T")
        # self.assertEqual(self.genome["i4000755"], "C")
        # self.assertEqual(self.genome["i4000759"], "G")

        self.assertEqual(self.genome["rs10488822"], "TC")
        self.assertEqual(self.genome["rs10810289"], "AA")
        self.assertEqual(self.genome["rs11980927"], "GG")
        self.assertEqual(self.genome["rs12913832"], "GG")
        self.assertEqual(self.genome["rs1426654"], "AA")
        self.assertEqual(self.genome["rs1540613"], "AG")
        self.assertEqual(self.genome["rs28504042"], "--")
        self.assertEqual(self.genome["rs3135027"], "G")
        self.assertEqual(self.genome["rs4477212"], "AT")
        self.assertEqual(self.genome["rs4536786"], "CA")
        self.assertEqual(self.genome["rs4672279"], "GT")
        self.assertEqual(self.genome["rs6015286"], "--")
        self.assertEqual(self.genome["rs6026400"], "CC")
        self.assertEqual(self.genome["rs6123756"], "TT")
        self.assertEqual(self.genome["rs742927"], "GG")
        self.assertEqual(self.genome["rs7715122"], "AT")
        self.assertEqual(self.genome["rs913897"], "AC")

        self.assertEqual(self.genome[10488822], "TC")
        self.assertEqual(self.genome[10810289], "AA")
        self.assertEqual(self.genome[11980927], "GG")
        self.assertEqual(self.genome[12913832], "GG")
        self.assertEqual(self.genome[1426654], "AA")
        self.assertEqual(self.genome[1540613], "AG")
        self.assertEqual(self.genome[28504042], "--")
        self.assertEqual(self.genome[3135027], "G")
        self.assertEqual(self.genome[4477212], "AT")
        self.assertEqual(self.genome[4536786], "CA")
        self.assertEqual(self.genome[4672279], "GT")
        self.assertEqual(self.genome[6015286], "--")
        self.assertEqual(self.genome[6026400], "CC")
        self.assertEqual(self.genome[6123756], "TT")
        self.assertEqual(self.genome[742927], "GG")
        self.assertEqual(self.genome[7715122], "AT")
        self.assertEqual(self.genome[913897], "AC")

        with self.assertRaises(TypeError):
            self.genome[1.0]

        with self.assertRaises(KeyError):
            self.genome[123]

        with self.assertRaises(KeyError):
            self.genome["rs123"]

        with self.assertRaises(KeyError):
            self.genome["rs91389"]

        with self.assertRaises(KeyError):
            self.genome["rs9138971"]

        with self.assertRaises(KeyError):
            self.genome["rs813897"]

        with self.assertRaises(KeyError):
            self.genome["rs13897"]

    def test_module_load_factor(self):
        self.assertIsInstance(self.genome.load_factor(), float)
        self.assertGreater(self.genome.load_factor(), 0.0)
        self.assertLess(self.genome.load_factor(), 1.0)
