import unittest

from starp.models import AmasPrimer, Sequence

class TestAmasPrimer(unittest.TestCase):

    def test_rev_comp(self):
        amas = AmasPrimer('GGCCGGTCCTAGTGTTGA', 1, (29, 47), strand=1)
        amas.tail = Sequence('GACGCAAGTGAGCAGTATGAC')

        revcomp = amas.rev_comp()

        expected = 'TCAACACTAGGACCGGCCGTCATACTGCTCACTTGCGTC'

        self.assertEqual(revcomp.strand, -1)
        self.assertEqual(str(revcomp), expected)

    def test_generate_amas_substitution(self):
        from starp.amasfactory import generate_amas_for_substitution
        # snp:                                      |
        allele1 = 'TGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAATGCTACGTGTGATCAACTTGGTCCACCACC'
        allele2 = 'TGTCGTTGTGGAGTCTATATCAACCACCATTCTAAAATGCTACGTGTGATCAACTTGGTCCACCACC'
        snp_position = 33

        upstream_amas, downstream_amas = generate_amas_for_substitution(
            allele1, allele2, snp_position
        )

        self.assertEqual('GAGTCTATATCAACCACCATTCTT', str(upstream_amas[0].sequence))
        self.assertEqual('GAGTCTATATCAACCACCATTCTA', str(upstream_amas[1].sequence))

        self.assertEqual('TAAATGCTACGTGTGATCAAC', str(downstream_amas[0].sequence))
        self.assertEqual('AAAATGCTACGTGTGATCAAC', str(downstream_amas[1].sequence))

if __name__ == "__main__":
    unittest.main()