"""
Tests for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import arv
import arv.match
import sys
import unittest

class ArvModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filename = "tests/fake_genome.txt"
        cls.genome = arv.load(cls.filename)

        cls.keys = sorted([
            "i3001754",
            "i3001755",
            "i3001759",
            "i3001761",
            "i3001773",
            "i4000755",
            "i4000759",
            "rs10488822",
            "rs10810289",
            "rs11980927",
            "rs12913832",
            "rs1426654",
            "rs1540613",
            "rs28504042",
            "rs3135027",
            "rs4477212",
            "rs4536786",
            "rs4672279",
            "rs6015286",
            "rs6026400",
            "rs6123756",
            "rs671",
            "rs742927",
            "rs7715122",
            "rs913897",
            ])

    def _int(self, rsid):
        """Converts RSID to integer."""
        if rsid.startswith("rs"):
            return int(rsid[2:])
        elif rsid.startswith("i"):
            return -int(rsid[1:])
        else:
            raise ValueError(rsid)

    def test_open_error(self):
        with self.assertRaises(RuntimeError):
            arv.load("non-existing-file")

    def test_len(self):
        self.assertEqual(len(self.genome), 25)

    def test_contains(self):
        for key in self.keys:
            self.assertIn(key, self.genome)
            self.assertNotIn("x" + key, self.genome)
            self.assertIn(self._int(key), self.genome)
            self.assertNotIn(-1*self._int(key), self.genome)

    def test_ychromosome(self):
        self.assertTrue(self.genome.y_chromosome)

    def test_ethnicity(self):
        self.assertEqual(self.genome.ethnicity, "")
        self.genome.ethnicity = "european"
        self.assertEqual(self.genome.ethnicity, "european")
        self.genome.ethnicity = ""
        self.assertEqual(self.genome.ethnicity, "")

    def test_orientation(self):
        self.assertEqual(self.genome.orientation, +1)
        self.genome.orientation = -1
        self.assertEqual(self.genome.orientation, -1)
        self.genome.orientation = +1
        self.assertEqual(self.genome.orientation, +1)
        with self.assertRaises(ValueError):
            self.genome.orientation = -2
        with self.assertRaises(ValueError):
            self.genome.orientation = 2

    def test_name(self):
        self.assertEqual(self.genome.name, self.filename)
        self.genome.name = "foo bar"
        self.assertEqual(self.genome.name, "foo bar")
        self.genome.name = self.filename
        self.assertEqual(self.genome.name, self.filename)

    def test_genotype_compare(self):
        a = self.genome["rs4477212"].genotype
        b = self.genome["rs4672279"].genotype

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

    def test_unphased_match(self):
        self.assertIsInstance(self.genome["rs10488822"], arv.SNP)
        self.assertEqual(self.genome["rs10488822"].genotype, "TC")

        self.assertEqual(arv.unphased_match(self.genome["rs10488822"], {
            "TC": "Matched TC",
            None: "No match"}), "Matched TC")

        self.assertEqual(arv.unphased_match(self.genome["rs10488822"], {
            "CT": "Matched CT",
            None: "No match"}), "Matched CT")

        self.assertEqual(arv.unphased_match(self.genome["rs10488822"], {
            "AT": "Matched AT",
            None: "No match"}), "No match")

        # No default case (None key)
        with self.assertRaises(KeyError):
            arv.unphased_match(self.genome["rs10488822"], {
                "AT": "Matched AT"})

    def test_snp_str_comparison(self):
        get = lambda key: self.genome[key]
        self.assertEqual(get("i3001754"), "A")
        self.assertEqual(get("i3001755"), "--")
        self.assertEqual(get("i3001759"), "--")
        self.assertEqual(get("i3001761"), "--")
        self.assertEqual(get("i3001773"), "T")
        self.assertEqual(get("i4000755"), "C")
        self.assertEqual(get("i4000759"), "G")
        self.assertEqual(get("rs10488822"), "TC")
        self.assertEqual(get("rs10810289"), "AA")
        self.assertEqual(get("rs11980927"), "GG")
        self.assertEqual(get("rs12913832"), "GG")
        self.assertEqual(get("rs1426654"), "AA")
        self.assertEqual(get("rs1540613"), "AG")
        self.assertEqual(get("rs28504042"), "--")
        self.assertEqual(get("rs3135027"), "G")
        self.assertEqual(get("rs4477212"), "AT")
        self.assertEqual(get("rs4536786"), "CA")
        self.assertEqual(get("rs4672279"), "GT")
        self.assertEqual(get("rs6015286"), "--")
        self.assertEqual(get("rs6026400"), "CC")
        self.assertEqual(get("rs6123756"), "TT")
        self.assertEqual(get("rs742927"), "GG")
        self.assertEqual(get("rs7715122"), "AT")
        self.assertEqual(get("rs913897"), "AC")

    def test_genotypes(self):
        get = lambda key: self.genome[key].genotype

        self.assertEqual(get("i3001754"), "A")
        self.assertEqual(get("i3001755"), "--")
        self.assertEqual(get("i3001759"), "--")
        self.assertEqual(get("i3001761"), "--")
        self.assertEqual(get("i3001773"), "T")
        self.assertEqual(get("i4000755"), "C")
        self.assertEqual(get("i4000759"), "G")
        self.assertEqual(get("rs10488822"), "TC")
        self.assertEqual(get("rs10810289"), "AA")
        self.assertEqual(get("rs11980927"), "GG")
        self.assertEqual(get("rs12913832"), "GG")
        self.assertEqual(get("rs1426654"), "AA")
        self.assertEqual(get("rs1540613"), "AG")
        self.assertEqual(get("rs28504042"), "--")
        self.assertEqual(get("rs3135027"), "G")
        self.assertEqual(get("rs4477212"), "AT")
        self.assertEqual(get("rs4536786"), "CA")
        self.assertEqual(get("rs4672279"), "GT")
        self.assertEqual(get("rs6015286"), "--")
        self.assertEqual(get("rs6026400"), "CC")
        self.assertEqual(get("rs6123756"), "TT")
        self.assertEqual(get("rs742927"), "GG")
        self.assertEqual(get("rs7715122"), "AT")
        self.assertEqual(get("rs913897"), "AC")

        self.assertEqual(get(-3001754), "A")
        self.assertEqual(get(-3001755), "--")
        self.assertEqual(get(-3001759), "--")
        self.assertEqual(get(-3001761), "--")
        self.assertEqual(get(-3001773), "T")
        self.assertEqual(get(-4000755), "C")
        self.assertEqual(get(-4000759), "G")
        self.assertEqual(get(10488822), "TC")
        self.assertEqual(get(10810289), "AA")
        self.assertEqual(get(11980927), "GG")
        self.assertEqual(get(12913832), "GG")
        self.assertEqual(get(1426654), "AA")
        self.assertEqual(get(1540613), "AG")
        self.assertEqual(get(28504042), "--")
        self.assertEqual(get(3135027), "G")
        self.assertEqual(get(4477212), "AT")
        self.assertEqual(get(4536786), "CA")
        self.assertEqual(get(4672279), "GT")
        self.assertEqual(get(6015286), "--")
        self.assertEqual(get(6026400), "CC")
        self.assertEqual(get(6123756), "TT")
        self.assertEqual(get(742927), "GG")
        self.assertEqual(get(7715122), "AT")
        self.assertEqual(get(913897), "AC")

        with self.assertRaises(KeyError):
            get(1.0)

        with self.assertRaises(KeyError):
            get(123)

        with self.assertRaises(KeyError):
            get(-123)

        with self.assertRaises(KeyError):
            get("rs123")

        with self.assertRaises(KeyError):
            get("rs123")

        with self.assertRaises(KeyError):
            get("i123")

        with self.assertRaises(KeyError):
            get("rs91389")

        with self.assertRaises(KeyError):
            get("rs9138971")

        with self.assertRaises(KeyError):
            get("rs813897")

        with self.assertRaises(KeyError):
            get("rs13897")

    def test_load_factor(self):
        self.assertIsInstance(self.genome.load_factor(), float)
        self.assertGreater(self.genome.load_factor(), 0.0)
        self.assertLess(self.genome.load_factor(), 1.0)

    def test_keys(self):
        # Depends on a static list
        self.assertEqual(sorted(list(self.genome.keys())), sorted(self.keys))

    def test_values(self):
        # Depends on genome.keys
        self.assertEqual(len(self.genome), len(list(self.genome.values())))
        self.assertEqual(len(self.keys), len(list(self.genome.values())))

        for rsid, snp in zip(self.genome.keys(), self.genome.values()):
            self.assertEqual(snp, self.genome[rsid])

    def test_items(self):
        # Depends on genome.keys and genome.values
        self.assertEqual(len(self.genome), len(list(self.genome.items())))
        self.assertEqual(len(self.keys), len(list(self.genome.items())))
        self.assertEqual(len(list(self.genome.keys())),
                len(list(self.genome.items())))
        self.assertEqual(len(list(self.genome.values())),
                len(list(self.genome.items())))

        for rsid, snp in self.genome.items():
            self.assertIn(rsid, self.keys)
            self.assertEqual(snp, self.genome[rsid])

    def test_snp_comparison(self):
        a = self.genome["rs4477212"]
        b = self.genome["rs4672279"]

        self.assertEqual(a, "AT")
        self.assertEqual(a.chromosome, 1)
        self.assertEqual(a.genotype, "AT")
        self.assertEqual(a.position, 82154)
        self.assertEqual(b, "GT")
        self.assertEqual(b.chromosome, 2)
        self.assertEqual(b.genotype, "GT")
        self.assertEqual(b.position, 59444675)
        self.assertGreater(b, a)
        self.assertGreaterEqual(a, a)
        self.assertGreaterEqual(b, a)
        self.assertGreaterEqual(b, b)
        self.assertLess(a, b)
        self.assertLessEqual(a, a)
        self.assertLessEqual(a, b)
        self.assertLessEqual(b, b)
        self.assertNotEqual(a, "GT")
        self.assertNotEqual(a, "TA")
        self.assertNotEqual(a, "TG")
        self.assertNotEqual(a, b)
        self.assertNotEqual(b, "")
        self.assertNotEqual(b, "AT")
        self.assertNotEqual(b, "TG")

    def test_snps(self):
        snp = self.genome.__getitem__

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
