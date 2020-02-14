import unittest
from starp.models import Snp

class TestSnp(unittest.TestCase):

    def test_simple_prefix(self):
        descriptor = '.12C>G'
        snp = Snp(descriptor)
        self.assertEqual(snp.prefix, '.')

    def test_init_one_subsitution(self):
        descriptor = 'NC_000023.10:g.33038255C>A'
        snp = Snp(descriptor)
        self.assertEqual(snp.prefix, 'NC_000023.10:g.')
        self.assertEqual(snp.position, 33038255)
        self.assertEqual(snp.type, 'substitution')
        self.assertEqual(snp.ref_nucleotide, 'C')
        self.assertEqual(snp.new_nucleotide, 'A')

    def test_init_one_insertion(self):
        descriptor = 'NC_000023.10:g.32867861insT'
        snp = Snp(descriptor)
        self.assertEqual(snp.prefix, 'NC_000023.10:g.')
        self.assertEqual(snp.position, 32867861)
        self.assertEqual(snp.type, 'insertion')
        self.assertEqual(snp.ref_nucleotide, '')
        self.assertEqual(snp.new_nucleotide, 'T')

    def test_init_one_deletion(self):
        descriptor = 'NG_012232.1:g.19del'
        snp = Snp(descriptor)
        self.assertEqual(snp.prefix, 'NG_012232.1:g.')
        self.assertEqual(snp.position, 19)
        self.assertEqual(snp.type, 'deletion')
        self.assertEqual(snp.ref_nucleotide, '')
        self.assertEqual(snp.new_nucleotide, '')

if __name__ == "__main__":
    unittest.main()