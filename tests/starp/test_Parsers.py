import unittest

from starp.models import Snp
from starp.parsers import TwoAlleles

class TestAmasPrimer(unittest.TestCase):

    def test_rev_comp(self):
        data = """>H26A region
TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTGCCTATTTTATTTA
>H26A region
TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGA-AATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTGCCTATTTTATTTA"""

        parser = TwoAlleles(data)
        snp = parser.snps()[0]

        expected = Snp(descriptor='.39del')
        expected.ref_nucleotide = 'A'
        
        self.assertEqual(expected, snp)
    

if __name__ == "__main__":
    unittest.main()