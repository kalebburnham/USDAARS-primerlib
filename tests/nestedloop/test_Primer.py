import unittest

from nestedloop.models import Sequence, Primer

class TestPrimer(unittest.TestCase):

    def test_self_complementary(self):
        seq = Sequence('AATT')
        primer = Primer(seq, 0, 0, 0)
        self.assertEqual(primer.complementary_score, 4)

    def test_self_complementary_2(self):
        """  
        AAAGGGGGT
          |     |
          TGGGGGAAA 
        """
        
        seq = Sequence('AAAGGGGGT')
        primer = Primer(seq, 0, 0, 0)
        self.assertEqual(primer.complementary_score, 2)

    def test_self_complementary_3(self):
        """  
           GGGATGGGGGG
              || 
        GGGGGGTAGGG 
        """
        seq = Sequence('GGGATGGGGGG')
        primer = Primer(seq, 0, 0, 0)
        self.assertEqual(primer.complementary_score, 2)

    def test_self_contig_complementary_score(self):
        """
        AAAACGCCTTTTT
        |||| |   ||||
        TTTTTCCGCAAAA
        """
        seq = Sequence('AAAACGCCTTTTT')
        primer = Primer(seq, 0, 0, 0)
        self.assertEqual(primer.contig_complementary_score, 4)

    def test_self_contig_complementary_score_2(self):
        """
        AATT
        ||||
        TTAA
        """
        seq = Sequence('AATT')
        primer = Primer(seq, 0, 0, 0)
        self.assertEqual(primer.complementary_score, 4)

if __name__ == "__main__":
    unittest.main()
