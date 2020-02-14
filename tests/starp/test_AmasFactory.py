import unittest

from starp.models import AmasPrimer, Sequence

class TestAmasPrimer(unittest.TestCase):

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
        from starp.amasfactory import generate_amas_for_indel

        allele1 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAATCGTCGT')
        allele2 = Sequence('TTGTCGTTGTGGAGTCTATATCAACCACCATTCTTA-TCGTCGT')
        snp_pos = 39

        upstream_amas, downstream_amas = generate_amas_for_indel(
                allele1, allele2, snp_pos
            )

if __name__ == "__main__":
    unittest.main()