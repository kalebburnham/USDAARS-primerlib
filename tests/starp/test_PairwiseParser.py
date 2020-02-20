import json
import unittest
from starp.parsers import PairwiseParser

class TestPairwiseParser(unittest.TestCase):

    data = """Score =  159 bits (86),  Expect = 1e-35
 Identities = 97/102 (95%), Gaps = 1/102 (1%)
 Strand=Plus/Minus

Query  1          TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  60
                  ||||||||  |||||||||||| |||||||||||||||||||||||||||||||||||||
Sbjct  528388742  TGCGATGA-TGAAAAAAAAAGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  528388684

Query  61         CCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGAT  102
                  ||| |||||||||||| |||||||||||||||||||||||||
Sbjct  528388683  CCGTGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  528388642


 Score =  156 bits (84),  Expect = 1e-34
 Identities = 95/100 (95%), Gaps = 2/100 (2%)
 Strand=Plus/Minus

Query  3          CGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACC  62
                  ||||||| | ||||||| || |||||||||||||||||||||||||||||||||||||||
Sbjct  147921921  CGATGAC-G-AAAAAAACGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACC  147921864

Query  63         GCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGAT  102
                  |||||||||||||| |||||||||||||||||||||||||
Sbjct  147921863  GCGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  147921824


 Score =  141 bits (76),  Expect = 3e-30
 Identities = 93/101 (92%), Gaps = 1/101 (1%)
 Strand=Plus/Minus

Query  1          TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  60
                  ||||||||| ||||||||| |  ||| |||||||||||||||||||||||||||||||||
Sbjct  144062963  TGCGATGAC-GAAAAAAAACGACGGTAGGAGTATGACGAAAATAAACCAGCGAAAATTAA  144062905

Query  61         CCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGA  101
                  |||||||||||||||| ||||||||||||||  ||||||||
Sbjct  144062904  CCGCGCGGACTATTCACCAAGTCGTCCATTAAAAGTAGAGA  144062864
"""
    def test_retrieved_hseqs(self):
        parser = PairwiseParser()
        hseqs = parser.parse(self.data)
        print(hseqs)
        expected = ['TGCGATGATGAAAAAAAAAGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGTGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT',
                    'CGATGACGAAAAAAACGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT',
                    'TGCGATGACGAAAAAAAACGACGGTAGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCACCAAGTCGTCCATTAAAAGTAGAGA']
        self.assertEqual(expected, hseqs)

if __name__ == "__main__":
    unittest.main()
