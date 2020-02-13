import unittest

from starp.models import AmasPrimer, Sequence

class TestAmasPrimer(unittest.TestCase):

    def test_rev_comp(self):
        amas = AmasPrimer('GGCCGGTCCTAGTGTTGA', 1, (29, 47), strand=1)
        amas.tail = Sequence('GACGCAAGTGAGCAGTATGAC')

        revcomp = amas.rev_comp()

        expected = 'TCAACACTAGGACCGGCCGTCATACTGCTCACTTGCGTC'
        print(amas.strand)
        print(revcomp.strand)
        self.assertEqual(revcomp.strand, -1)
        self.assertEqual(str(revcomp), expected)

if __name__ == "__main__":
    unittest.main()