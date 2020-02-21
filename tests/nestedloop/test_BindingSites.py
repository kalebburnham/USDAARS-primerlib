import unittest

from nestedloop.models import Sequence, Primer
from nestedloop.utils import binding_sites

class TestBindingSites(unittest.TestCase):

    def test_no_binding_sites(self):
        ref = Sequence('GTAGCATGCGATCGTCTATCGTATGACGCTGCACGCAT')
        primer = Primer('TGACGGAAAAAAA', span=(13, 26), strand=1)
        sites = binding_sites([ref], primer)
        expected = []
        self.assertEqual(expected, [site.span() for site in sites])

    def test_single_binding_site(self):
        ref = Sequence('TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTGCCTATTTTATTTATTCATTGGAAAACTTGAAGAATGGAAATTTTGTGGTGCGATAATAGCATATTGTGTTGATGTTGTCAAATTTCTTCCCACCAACCGTGCCCTCCTCCAGTTGATTGTATGCACAAGCGTGCCACCCCATGACGGTAGGAAACGAGCACAACTGCACCATTGTTTAAATCATCTCGTCAGAGCTTCATCAGTTTATAATGGCATATGTTTGTCCACAAGCGCATCACCGGCCCCGGTGCGCCAAATATATTTTGTAACCACCCTTATCCTTTCAATCAACCATCACGAGGAAATAGTCTGAGAAGTTCGCAGCCACATTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAAATTAGCTACGTGTGATCAACTTGGTCCACCACCTACGCAAGGAGGTGCGTGGTGCGGAAGTTGATTAATCCAACCAAACGTGCCTGCTTGTCTCTAGATAATACAAACTAAGATGTCTAACTAAACTAATGCATCTGCTCAAAGGCATCAATCGCATCGGCGCAAGCGTTCAGAGGTCAGACGACCATTGGCAAGCACCAAAGGACGGGGCCATACGCTTAGAGACACCTACATGGCCGGGTATTAATTACTAGAGTGTTGCAACGCGTTCGCAAAGTCAGCACACGGAGATTGGTGC')
        primer = Primer('CAGTTGATTGTATGCACAAGCGTGCCACCCCATGA', span=(228, 263), strand=1)
        sites = binding_sites([ref], primer)
        expected = [(228, 263)]
        self.assertEqual(expected, [site.span() for site in sites])

    def test_multiple_binding_sites(self):
        ref = Sequence('TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGACGGAAAAAAAAAGGTGGTGGGAGTATGACGGATTTTCACAACCCAATTTGCCTATTTTATTTATTCATTGGAAAACTTGAAGAATGGAAATTTTGTGGTGCGATAATAGCATATTGTGTTGATGTTGTC')
        primer = Primer('GACGGAAAAAAAAAGGTGGTGGGAGTATGACG', span=(6, 38), strand=1)
        sites = binding_sites([ref], primer)
        expected = [(6, 38), (99, 131)]
        self.assertEqual(expected, [site.span() for site in sites])

    def test_fuzzy_binding_sites(self):
        # ACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAAT
        # ||||  |||||||| |||||| ||||||||||||||
        # ACGGTTAAAAAAAGATGGTGGCAGTATGACGAAAAT
        ref = Sequence('TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGACGGAAAAAAAAAGGTGGTGGGAGTATGACGGACGGTTAAAAAAAGATGGTGGCAGTATGACGAAAATATTTTCACAACCCAATTTGCCTATTTTATTTATTCATTGGAAAACTTGAAGAATGGAAATTTTGTGGTGCGATAATAGCATATTGTGTTGATGTTGTC')
        primer = Primer('ACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAAT', span=(6, 38), strand=1)
        sites = binding_sites([ref], primer)
        expected = [(7, 43), (132, 168)]
        self.assertEqual(expected, [site.span() for site in sites])

    def test_rev_comp_binding_sites(self):
        # TGACGGAAAAAAAAAGGTGGTGGGAGTA
        # Reverse Complement:
        # TACTCCCACCACCTTTTTTTTTCCGTCA
        ref = Sequence('TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTTACTCCCACCACCTTTTTTTTTCCGTCACATCAAGTCGTCCATTAGGAGTAG')
        primer = Primer('TGACGGAAAAAAAAAGGTGGTGGGAGTA', span=(5, 33), strand=1)
        sites = binding_sites([ref], primer)
        expected = [(5, 33), (74, 102)]
        self.assertEqual(expected, [site.span() for site in sites])

    def test_multiple_references(self):
        ref1 = Sequence('TGCGATGACGCGTGAGGTGTATGAAGCAAAGGCAGCGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGCGCGGACTATTTACTCCCACCACCTTTTTTTTTCCGTCACATCAAGTCGTCCATTAGGAGTAG')
        ref2 = Sequence('TATAAGATTATGAACGTGAGGTGTATGAAGCAAAGGCAGCAATCAAGCTAATATGGAGATGGATTTCAGAAGGTTTTGTTTATAGCGAAAAACAAGAAACTGACTTATATGAGCTCGGTAAGAAGTACTTCAATGAGCTAGTAAATAGAAGTATGATACAGCCAATTGGTATTGATGATGGAGAAGATAAACAAGCGTGTCGTGTACATGACATGGTGCTTGATATCCTATGCT')
        primer = Primer('CGTGAGGTGTATGAAGCAAAGGCAGC', span=(10, 36), strand=1)
        sites = binding_sites([ref1, ref2], primer)
        expected = [(10, 36), (14, 40)]
        self.assertEqual(expected, [site.span() for site in sites])

if __name__ == "__main__":
    unittest.main()
