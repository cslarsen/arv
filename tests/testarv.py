import arv
import os
import unittest

class ArvTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = arv.load("fake_genome.txt")

    def test_len(self):
        self.assertEqual(len(self.genome), 12)

    def test_ychromo(self):
        self.assertTrue(self.genome.ychromo())
