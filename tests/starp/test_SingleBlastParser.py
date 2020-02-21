import unittest

from starp.models import Snp
from starp.parsers import SingleBlastParser

class TestSingleBlastParser(unittest.TestCase):

    def test_single_blast_parser(self):
        sequence = """
Query  1          TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  60
                  |||||||| ||||||||||||| |||||||||||||||||||||||||||||||||||||
Sbjct  528388742  TGCGATGA-GGAAAAAAAAAGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  528388684

Query  61         CCG-GCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  102
                  ||| ||||||||||||||||||||||||||||||||||||||
Sbjct  528388683  CCGTGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  528388642"""
        parser = SingleBlastParser(sequence)
        breakpoint()
        snps = parser.snps()
        expected = [Snp('.8del'), Snp('.22T>C'), Snp('.63insT')]
        expected[0].ref_nucleotide = 'C'
        self.assertEqual(expected, snps)

if __name__ == "__main__":
    unittest.main()
