"""
Inferring tests for arv.

arv
Copyright 2017 Christian Stigen Larsen
Distributed under the GNU GPL v3 or later; see COPYING.
"""

import arv
import unittest

class ArvInferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("tests/fake_genome.txt")

    def test_infer_gender(self):
        gender = "man" if self.genome.y_chromosome else "woman"
        self.assertEqual(gender, "man")

    def test_infer_complexion(self):
        complexion = "light" if self.genome["rs1426654"] == "AA" else "dark"
        self.assertEqual(complexion, "light")

    def test_infer_unphased_match_eyecolor(self):
        eyecolor = arv.unphased_match(self.genome["rs12913832"], {
            "AA": "brown eyes",
            "AG": "brown or green eyes",
            "GG": "blue eyes"})
        self.assertEqual(eyecolor, "blue eyes")
