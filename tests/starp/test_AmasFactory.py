import unittest

from starp.models import AmasPrimer, Sequence, Snp

class TestAmasFactory(unittest.TestCase):

    def test_generate_amas_downstream_sub(self):
        """ Verify the spans are the same after generating downstream. """
        from starp.amasfactory import generate_amas_downstream

        allele1 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAT')
        allele2 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAT')

        smallest_primer1 = generate_amas_downstream(allele1, 1, 30, 5, 10)[0]
        smallest_primer2 = generate_amas_downstream(allele2, 2, 30, 5, 10)[0]

        expected1 = AmasPrimer(sequence='TTCTT', allele_num=1, span=(30, 35), strand=1)
        expected2 = AmasPrimer(sequence='TTCTT', allele_num=2, span=(30, 35), strand=1)

        self.assertEqual(expected1.span, smallest_primer1.span)
        self.assertEqual(expected2.span, smallest_primer2.span)

    def test_generate_amas_downstream_del(self):
        """ Indels require some manipulation to keep the spans correct.
        """
        from starp.amasfactory import generate_amas_downstream

        allele1 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAT')
        allele2 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTA-T')

        smallest_primer1 = generate_amas_downstream(allele1, 1, 30, 5, 10)[0]
        smallest_primer2 = generate_amas_downstream(allele2, 2, 30, 5, 10)[0]

        expected1 = AmasPrimer(sequence='TTCTT', allele_num=1, span=(30, 35), strand=1)
        expected2 = AmasPrimer(sequence='TTCTT', allele_num=2, span=(30, 35), strand=1)

        self.assertEqual(expected1.span, smallest_primer1.span)
        self.assertEqual(expected2.span, smallest_primer2.span)

    def test_generate_amas_for_indel(self):
        # TODO
        from starp.amasfactory import generate_amas_for_indel

        allele1 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAATCGTCGT')
        allele2 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTA-TCGTCGT')
        snp_pos = 39

        upstream_amas, downstream_amas = generate_amas_for_indel(
                allele1, allele2, snp_pos
        )

    def test_substitute_bases_one_snp_first(self):
        """ Since the snp position is at the first nucleotide, the
        sequences need to be reverse complemented THEN have its
        bases substituted. Order is important here. 
        
        These substitution tests are not conclusive. There are way too
        many branches to test every possible combination of SNPs.
        Rather, these tests show that the logic is at least somewhat
        functional. """
        from starp.amasfactory import substitute_bases

        seq1 = Sequence('TGCCTAACATAGAAAACATTGAGGCA')
        seq2 = Sequence('AGCCTAACATAGAAAACATTGAGGCA')
        snp_position = 'first'

        expected = (Sequence('TGCTTAACATAGAAAACATTGAGGCA'),
                    Sequence('AGTCTAACATAGAAAACATTGAGGCA'))

        self.assertEqual(expected, substitute_bases((seq1, seq2), snp_position))

    def test_substitute_bases_one_snp_last(self):
        from starp.amasfactory import substitute_bases

        seq1 = Sequence('ACGGAGTTACAAAAGATACAATCCGT')
        seq2 = Sequence('ACGGAGTTACAAAAGATACAATCCGA')
        snp_position = 'last'

        expected = (Sequence('ACGGAGTTACAAAAGATACAATCTGT'),
                    Sequence('ACGGAGTTACAAAAGATACAATTCGA'))

        self.assertEqual(expected, substitute_bases((seq1, seq2), snp_position))

    def test_substitute_bases_two_snps_first(self):
        """ Since the snp position is at the first nucleotide, the
        sequences need to be reverse complemented THEN have its
        bases substituted. Order is important here. """
        from starp.amasfactory import substitute_bases

        seq1 = Sequence('TGCCTAACATAGAAAACATTGAGGCA')
        seq2 = Sequence('ACCCTAACATAGAAAACATTGAGGCA')
        snp_position = 'first'

        expected = (Sequence('TGCTTAACATAGAAAACATTGAGGCA'),
                    Sequence('ACTCTAACATAGAAAACATTGAGGCA'))

        self.assertEqual(expected, substitute_bases((seq1, seq2), snp_position))

    def test_substitute_bases_two_snps_last(self):
        from starp.amasfactory import substitute_bases

        seq1 = Sequence('ACGGAGTTACAAAAGATACAATCCGT')
        seq2 = Sequence('ACGGAGTTACAAAAGATACAATCCTA')
        snp_position = 'last'

        expected = (Sequence('ACGGAGTTACAAAAGATACAATCTGT'),
                    Sequence('ACGGAGTTACAAAAGATACAATTCTA'))

        self.assertEqual(expected, substitute_bases((seq1, seq2), snp_position))

if __name__ == "__main__":
    unittest.main()