import unittest

from starp.models import Snp, Sequence
from starp.parsers import get_parser, TwoAlleles, SnpSequence
from starp.exceptions import StarpError

class TestFormatType(unittest.TestCase):
    def test_SnpSequence(self):
        data = 'TGCACAGTGTACTAGC[G/A]GTCGTATGACTCA'
        parser = get_parser(data)
        self.assertEqual(type(parser), SnpSequence)

    def test_TwoAlleles(self):
        data = ">Allele1\nTGACACGTACGT\n>Allele2\nTGACACGTACGT"
        parser = get_parser(data)
        self.assertEqual(type(parser), TwoAlleles)

class TestSnpSequence(unittest.TestCase):
    def test_NoSnps(self):
        data = "TGACACGTACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        self.assertEqual(len(snps), 0)

    def test_Substitution(self):
        # Substitution at position 5.
        data = "TGACA[C/G]GTACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.5C>G')]
        self.assertEqual(snps, expected)

    def test_Two_Substitutions(self):
        # Substitution at position 5 and 6.
        data = "TGACA[C/G][G/A]TACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.5C>G'), Snp('.6G>A')]
        self.assertEqual(snps, expected)

    def test_One_Insertion(self):
        # Insertion of 'A' at position 1.
        data = "T[-/A]GACACGTACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.1insA')]
        self.assertEqual(snps, expected)

    def test_Two_Insertions(self):
        # Allele1: TGACACGTACGT
        # Insertion of 'A' at position 1 and
        # Insertion of 'G' at position 6.
        data = "T[-/A]GACA[-/G]CGTACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.1insA'), Snp('.6insG')]
        self.assertEqual(snps, expected)

    def test_Two_Insertions_Consecutive(self):
        # Two insertions at position 1 and 2.
        data = "T[-/A][-/G]GACACGTACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.1insA'), Snp('.2insG')]
        self.assertEqual(snps, expected)

    def test_Two_Deletions(self):
        # Two deletions at positions 0 and 6.
        data = "[T/-]GACAC[G/-]TACGT"
        parser = SnpSequence(data)
        snps = parser.snps()
        expected = [Snp('.0del'), Snp('.6del')]
        expected[0].ref_nucleotide = 'T'
        expected[1].ref_nucleotide = 'G'
        self.assertEqual(snps, expected)

    def test_SyntaxError(self):
        # Raise Starp Error if a bad string is entered.
        data = "TGAC[ACGTACGT"
        with self.assertRaises(StarpError):
            SnpSequence(data)

    def test_InvalidChars(self):
        # Raise Starp Error if there are bad characters.
        data = "TGAC999ACGTACGT"
        with self.assertRaises(StarpError):
            SnpSequence(data)

    def test_SnpSameCharacter(self):
        # Raise Starp Error if Snps have the same character.
        data = "TGAC[A/A]CGTACGT"
        with self.assertRaises(StarpError):
            SnpSequence(data)

class TestTwoAlleles(unittest.TestCase):

    def test_NoSnps(self):
        data = ">Allele1\nTGACACGTACGT\n>Allele2\nTGACACGTACGT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        self.assertEqual(len(snps), 0)

    def test_Substitution(self):
        # Substitution at position 6: G>C
        data = ">Allele1\nTGACACGTACGT\n>Allele2\nTGACACCTACGT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.6G>C')]
        self.assertEqual(snps, expected)

    def test_One_Insertion(self):
        data = ">Allele1\nT-GACACGTACGT\n>Allele2\nTAGACACGTACGT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.1insA')]
        self.assertEqual(snps, expected)

    def test_Two_Insertions(self):
        # Allele1 = TGACACGTACGT
        # Insertion of 'A' at position 1 and 'G' at position 8.
        data = ">Allele1\nT-GACACG-TACGT\n>Allele2\nTAGACACGGTACGT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.1insA'), Snp('.8insG')]
        self.assertEqual(snps, expected)

    def test_One_Deletion(self):
        # Deletion at position 9.
        data = ">Allele1\nTGACACGTACGT\n>Allele2\nTGACACGTA-GT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.9del')]
        expected[0].ref_nucleotide = 'C'
        self.assertEqual(snps, expected)

    def test_Two_Deletions(self):
        # Two deletions at positions 9 and 11.
        data = ">Allele1\nTGACACGTACGT\n>Allele2\nTGACACGTA-G-"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.9del'), Snp('.11del')]
        expected[0].ref_nucleotide = 'C'
        expected[1].ref_nucleotide = 'T'
        self.assertEqual(snps, expected)

    def test_del_and_ins(self):
        # Deletion at position 0 and insertion of 'T' at position 4.
        data = ">Allele1\nTGAC-ACGTACGT\n>Allele2\n-GACTACGTACGT"
        parser = TwoAlleles(data)
        snps = parser.snps()
        expected = [Snp('.0del'), Snp('.4insT')]
        expected[0].ref_nucleotide = 'T'
        self.assertEqual(snps, expected)

    def test_invalid_chars(self):
        # StarpError should be thrown when there are invalid chars
        data = ">Allele1\nTGACAC5TACGT\n>Allele2\nTGACACGTACGT"
        with self.assertRaises(StarpError):
            TwoAlleles(data)

    def test_unequal_length(self):
        # StarpError should be thrown when alleles have unequal length.
        data = ">Allele1\nTGACACGTA\n>Allele2\nTGACACGTACGT"
        with self.assertRaises(StarpError):
            TwoAlleles(data)

class TestAmasPrimer(unittest.TestCase):

    def test_rev_comp(self):
        data = """>H26A region
TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTGCCTATTTTATTTA
>H26A region
TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGA-AATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTGCCTATTTTATTTA"""

        parser = TwoAlleles(data)
        snp = parser.snps()[0]

        expected = Snp(descriptor='.39del')
        expected.ref_nucleotide = 'A'
        
        self.assertEqual(expected, snp)
    

if __name__ == "__main__":
    unittest.main()