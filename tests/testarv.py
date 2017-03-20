"""
Tests for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import arv
import unittest

class ArvModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("tests/fake_genome.txt")

    def test_len(self):
        self.assertEqual(len(self.genome), 24)
        self.assertEqual(len(arv.Genome()), 0)

    def test_ychromosome(self):
        self.assertTrue(self.genome.y_chromosome)
        self.assertFalse(arv.Genome().y_chromosome)

    def test_genotype_compare(self):
        a = self.genome.get_snp("rs4477212").genotype
        b = self.genome.get_snp("rs4672279").genotype

        self.assertFalse(a != a)
        self.assertTrue(a != b)
        self.assertTrue(a < b)
        self.assertTrue(a <= b)
        self.assertTrue(a == a)
        self.assertTrue(b != a)
        self.assertTrue(b == b)
        self.assertTrue(b > a)
        self.assertTrue(b >= a)

        self.assertFalse("AT" > a)
        self.assertFalse(a > "AT")
        self.assertTrue("AT" >= a)
        self.assertTrue(a != "AA")
        self.assertTrue(a != "TA")
        self.assertTrue(a < b)
        self.assertTrue(a <= b)
        self.assertTrue(a == "AT")
        self.assertTrue(a > "AA")
        self.assertTrue(a > "AS")
        self.assertTrue(a >= "AA")
        self.assertTrue(a >= "AS")
        self.assertTrue(a >= "AT")
        self.assertTrue(b == "GT")
        self.assertTrue(b > "AT")
        self.assertTrue(b >= "AT")

    def test_genotypes(self):
        self.assertEqual(self.genome["i3001754"], "A")
        self.assertEqual(self.genome["i3001755"], "--")
        self.assertEqual(self.genome["i3001759"], "--")
        self.assertEqual(self.genome["i3001761"], "--")
        self.assertEqual(self.genome["i3001773"], "T")
        self.assertEqual(self.genome["i4000755"], "C")
        self.assertEqual(self.genome["i4000759"], "G")
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

        self.assertEqual(self.genome[-3001754], "A")
        self.assertEqual(self.genome[-3001755], "--")
        self.assertEqual(self.genome[-3001759], "--")
        self.assertEqual(self.genome[-3001761], "--")
        self.assertEqual(self.genome[-3001773], "T")
        self.assertEqual(self.genome[-4000755], "C")
        self.assertEqual(self.genome[-4000759], "G")
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

        with self.assertRaises(KeyError):
            self.genome[1.0]

        with self.assertRaises(KeyError):
            self.genome[123]

        with self.assertRaises(KeyError):
            self.genome[-123]

        with self.assertRaises(KeyError):
            self.genome["rs123"]

        with self.assertRaises(KeyError):
            self.genome["i123"]

        with self.assertRaises(KeyError):
            self.genome["rs91389"]

        with self.assertRaises(KeyError):
            self.genome["rs9138971"]

        with self.assertRaises(KeyError):
            self.genome["rs813897"]

        with self.assertRaises(KeyError):
            self.genome["rs13897"]

    def test_load_factor(self):
        self.assertIsInstance(self.genome.load_factor(), float)
        self.assertGreater(self.genome.load_factor(), 0.0)
        self.assertLess(self.genome.load_factor(), 1.0)

    def test_snps(self):
        snp = self.genome.get_snp

        self.assertEqual(snp("rs4477212").chromosome, 1)
        self.assertEqual(snp("rs4477212").position, 82154)
        self.assertEqual(snp("rs4477212").genotype, "AT")

        self.assertEqual(snp("rs4672279").chromosome, 2)
        self.assertEqual(snp("rs4672279").position, 59444675)
        self.assertEqual(snp("rs4672279").genotype, "GT")

        self.assertEqual(snp("rs4536786").chromosome, 3)
        self.assertEqual(snp("rs4536786").position, 140049121)
        self.assertEqual(snp("rs4536786").genotype, "CA")

        self.assertEqual(snp("rs7715122").chromosome, 5)
        self.assertEqual(snp("rs7715122").position, 94197884)
        self.assertEqual(snp("rs7715122").genotype, "AT")

        self.assertEqual(snp("rs11980927").chromosome, 7)
        self.assertEqual(snp("rs11980927").position, 20010422)
        self.assertEqual(snp("rs11980927").genotype, "GG")

        self.assertEqual(snp("rs10810289").chromosome, 9)
        self.assertEqual(snp("rs10810289").position, 14899708)
        self.assertEqual(snp("rs10810289").genotype, "AA")

        self.assertEqual(snp("rs10488822").chromosome, 11)
        self.assertEqual(snp("rs10488822").position, 35984271)
        self.assertEqual(snp("rs10488822").genotype, "TC")

        self.assertEqual(snp("rs913897").chromosome, 13)
        self.assertEqual(snp("rs913897").position, 73892459)
        self.assertEqual(snp("rs913897").genotype, "AC")

        self.assertEqual(snp("rs1540613").chromosome, 16)
        self.assertEqual(snp("rs1540613").position, 80476182)
        self.assertEqual(snp("rs1540613").genotype, "AG")

        self.assertEqual(snp("rs6123756").chromosome, 20)
        self.assertEqual(snp("rs6123756").position, 56556146)
        self.assertEqual(snp("rs6123756").genotype, "TT")

        self.assertEqual(snp("rs6015286").chromosome, 20)
        self.assertEqual(snp("rs6015286").position, 57048415)
        self.assertEqual(snp("rs6015286").genotype, "--")

        self.assertEqual(snp("rs6026400").chromosome, 20)
        self.assertEqual(snp("rs6026400").position, 57183524)
        self.assertEqual(snp("rs6026400").genotype, "CC")

        self.assertEqual(snp("rs742927").chromosome, "Y")
        self.assertEqual(snp("rs742927").position, 57183914)
        self.assertEqual(snp("rs742927").genotype, "GG")

        self.assertEqual(snp("i3001754").chromosome, "MT")
        self.assertEqual(snp("i3001754").position, 16256)
        self.assertEqual(snp("i3001754").genotype, "A")

        self.assertEqual(snp("i3001755").chromosome, "MT")
        self.assertEqual(snp("i3001755").position, 16257)
        self.assertEqual(snp("i3001755").genotype, "--")

        self.assertEqual(snp("i3001759").chromosome, "MT")
        self.assertEqual(snp("i3001759").position, 16258)
        self.assertEqual(snp("i3001759").genotype, "--")

        self.assertEqual(snp("i3001761").chromosome, "MT")
        self.assertEqual(snp("i3001761").position, 16259)
        self.assertEqual(snp("i3001761").genotype, "--")

        self.assertEqual(snp("i3001773").chromosome, "MT")
        self.assertEqual(snp("i3001773").position, 16265)
        self.assertEqual(snp("i3001773").genotype, "T")

        self.assertEqual(snp("i4000755").chromosome, "MT")
        self.assertEqual(snp("i4000755").position, 16548)
        self.assertEqual(snp("i4000755").genotype, "C")

        self.assertEqual(snp("i4000759").chromosome, "MT")
        self.assertEqual(snp("i4000759").position, 16567)
        self.assertEqual(snp("i4000759").genotype, "G")

        self.assertEqual(snp("rs1426654").chromosome, 15)
        self.assertEqual(snp("rs1426654").position, 48426484)
        self.assertEqual(snp("rs1426654").genotype, "AA")

        self.assertEqual(snp("rs12913832").chromosome, 15)
        self.assertEqual(snp("rs12913832").position, 28365618)
        self.assertEqual(snp("rs12913832").genotype, "GG")

        self.assertEqual(snp("rs28504042").chromosome, "MT")
        self.assertEqual(snp("rs28504042").position, 1549)
        self.assertEqual(snp("rs28504042").genotype, "--")

        self.assertEqual(snp("rs3135027").chromosome, "MT")
        self.assertEqual(snp("rs3135027").position, 1598)
        self.assertEqual(snp("rs3135027").genotype, "G")
