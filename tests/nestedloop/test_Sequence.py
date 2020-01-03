import unittest

from nestedloop.models import Sequence

class TestSequenceMethods(unittest.TestCase):

    def test_complement(self):
        seq = Sequence('ACGT')
        expected = Sequence('TGCA')
        self.assertEqual(seq.complement(), expected)

    def test_complement_invalid_chars(self):
        seq = Sequence('ACGTN')
        expected = Sequence('TGCAN')
        self.assertEqual(seq.complement(), expected)

    def test_gc(self):
        seq = Sequence('GAGATCTC')
        expected = 0.50
        self.assertEqual(seq.gc, expected)

    def test_gc_empty(self):
        """ gc() divides by 0 """
        seq = Sequence('')
        expected = 0
        self.assertEqual(seq.gc, expected)

if __name__ == "__main__":
    unittest.main()