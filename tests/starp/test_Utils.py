import unittest

from starp.models import Primer, Sequence
from starp.utils import cut, rsorted

class TestUtils(unittest.TestCase):

    def test_two_primer_sort(self):
        good_primer = Primer('TGCGCTTTCTC', None, None, 1)
        bad_primer = Primer('AAAAAAAAAAA', None, None, 1)
        primers = [good_primer, bad_primer]
        expected = [good_primer, bad_primer]

        self.assertEqual(expected, rsorted(primers))

    def test_cut(self):
        tail = Sequence('GCAACAGGAACCAGCTATGAC')
        primer = Primer('CTCACACCTCTCAAACGACTG', (0, 0), (0, 0), 1)

        expected = Sequence('GCAACAGGAACCAGCTATGA')
        self.assertEqual(expected, cut(tail, primer))

if __name__ == "__main__":
    unittest.main()