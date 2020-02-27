import unittest

from starp.models import Primer
from starp.utils import rsorted

class TestUtils(unittest.TestCase):

    def test_two_primer_sort(self):
        good_primer = Primer('TGCGCTTTCTC', None, None, 1)
        bad_primer = Primer('AAAAAAAAAAA', None, None, 1)
        primers = [good_primer, bad_primer]
        expected = [good_primer, bad_primer]

        self.assertEqual(expected, rsorted(primers))

if __name__ == "__main__":
    unittest.main()