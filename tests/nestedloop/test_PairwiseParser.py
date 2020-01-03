import sys
sys.path.append('../..')

import json
import unittest
from NestedLoop.Parsers import PairwiseParser

class TestPairwiseParser(unittest.TestCase):

    ncbi_test_results = """Triticum monococcum subsp. monococcum cultivar DV92 Sr35 region, genomic sequence 
Sequence ID: KC573058.1Length: 307519Number of Matches: 8
Range 1: 18021 to 18405GenBankGraphics
Next Match
Previous Match
Alignment statistics for match #1
Score
Expect
Identities
Gaps
Strand
145 bits(160)
3e-30
270/391(69%)
48/391(12%)
Plus/Plus
Query  433    GAAGTTCGCAGCCACATTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAAATTAGC  492
              |||||| || ||   | | | | || |||||||| |||||||||| ||||||||||||||
Sbjct  18021  GAAGTTAGCTGCGGGAATATAGCTGCGGAGTCTAGATCAACCACCGTTCTTAAAATTAGC  18080

Query  493    TACGTGTGATCAACTTGGTCCACCACCTACGCAAGGAGGTGCGTGGTGCGGAAGTTGATT  552
               |  |||| || | ||| ||||| ||  ||||||||||||||||||| |  |||||||||
Sbjct  18081  AATATGTGGTC-ATTTGCTCCACTACGGACGCAAGGAGGTGCGTGGTTCCAAAGTTGATT  18139

Query  553    -AATCCAAC--CAAACGTGCCTGCTTGTCTCTAGATAATAC-AAACTAAGATGTCTAACT  608
               ||||||||  | || |   |||   || |||||| ||| | ||||| || | | |||||
Sbjct  18140  AAATCCAACGACCAAGGCAGCTGTGGGTTTCTAGACAATGCAAAACTTAGCTATTTAACT  18199

Query  609    AAACTAATGCATCTGCTCA---------------AAGGCATCAATCGCATCGGCGCAAGC  653
              | ||||||||| ||| | |                || |||| || |||||  | |||| 
Sbjct  18200  ACACTAATGCACCTGTTTAGAGAGTTGGGGCATCGAGCCATCGATGGCATCATCTCAAGT  18259

Query  654    GTTCAGAGGTCAGACGACCATTGGCAAGCACCAAAGGACGGGGC--------CATACGCT  705
              || | || | |    | | |||||||||||||||||||| ||||        ||||||||
Sbjct  18260  GTCCTGACGCCCATTGTCTATTGGCAAGCACCAAAGGACAGGGCCAAAGGGGCATACGCT  18319

Query  706    TA---GAG-ACACCTACATGGCCGGGTATTAATTA-------CTAGA-GTGTTGCAACGC  753
              ||    ||  || | |||||||    |||||||||       | | | | ||  ||||||
Sbjct  18320  TATCCTAGCGCATCCACATGGC----TATTAATTAGTGTGTTCCACATGGGTGACAACGC  18375

Query  754    GTTCGCAAAGTCAGCACACGGA---GATTGG  781
              | ||||||||||||| || |||   ||||||
Sbjct  18376  G-TCGCAAAGTCAGCGCATGGAATCGATTGG  18405


Range 2: 17154 to 17429GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #2
Score
Expect
Identities
Gaps
Strand
136 bits(150)
2e-27
205/287(71%)
21/287(7%)
Plus/Plus
Query  121    CTATTTTATTTATTCATTGGAAAACTTGAA----GAATGGAAATTTTGTGGTGCGATAAT  176
              |||||||||||||| ||  |||| ||||||    || | |||| | ||||||| ||||| 
Sbjct  17154  CTATTTTATTTATTGATGTGAAA-CTTGAATCTAGATTAGAAACTGTGTGGTGTGATAAA  17212

Query  177    AGCATATTGTGTTGATGTTGTCAAATTTCTTCCCACCAACCGTGCCCTCCTCCA-----G  231
              ||||||| |||||||||||||||||||||| |||       |||  | || | |     |
Sbjct  17213  AGCATATCGTGTTGATGTTGTCAAATTTCTGCCCTTTTTT-GTGAACACCACAAACCGTG  17271

Query  232    TTGATTGTA-TGCACAAGCGTGCCACCCCATGACGGTAGGAAACGAGCACAACTGCACCA  290
              |||||  || | ||      | |||||   ||| ||  |   | ||||| ||   ||  |
Sbjct  17272  TTGATGATACTTCA------TCCCACCAACTGATGGCTG---ATGAGCAGAATGACAAAA  17322

Query  291    TTGTTTAAATCATCTCGTCAGAGCTTCATCAGTTTATAATGGCATATGTTTGTCCACAAG  350
              |||||| |||  | || | ||  |||| ||||||||| ||||||| ||||||||||||||
Sbjct  17323  TTGTTTTAATTCTATCATGAGGACTTCTTCAGTTTATGATGGCATCTGTTTGTCCACAAG  17382

Query  351    CGCATCACCGGCCCCGGTGCGCCAAATATATTTTGTAACCACCCTTA  397
              ||||| |   || | | |||  ||||||||||||||| ||| |||||
Sbjct  17383  CGCATTATATGCTCTGATGCAGCAAATATATTTTGTATCCATCCTTA  17429


Range 3: 61107 to 61216GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #3
Score
Expect
Identities
Gaps
Strand
68.0 bits(74)
7e-07
83/111(75%)
2/111(1%)
Plus/Plus
Query  546    GTTGATTAATCCAACCAAACGTGCCTGCTTGTCTCTAGATAATACAAACTAAGATGTCTA  605
              |||||||| |||||||||   ||  ||| |  ||||||||||| |||||  || |||  |
Sbjct  61107  GTTGATTAGTCCAACCAAGGTTGT-TGCCTACCTCTAGATAATGCAAACCTAGCTGTTAA  61165

Query  606    ACTA-AACTAATGCATCTGCTCAAAGGCATCAATCGCATCGGCGCAAGCGT  655
              ||||   ||| ||| |||||||| ||| | | | | ||||  |||||||||
Sbjct  61166  ACTAGCGCTAGTGCGTCTGCTCAGAGGTAACGAGCACATCCTCGCAAGCGT  61216


Range 4: 40301 to 40332GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #4
Score
Expect
Identities
Gaps
Strand
59.0 bits(64)
3e-04
32/32(100%)
0/32(0%)
Plus/Plus
Query  604    TAACTAAACTAATGCATCTGCTCAAAGGCATC  635
              ||||||||||||||||||||||||||||||||
Sbjct  40301  TAACTAAACTAATGCATCTGCTCAAAGGCATC  40332


Range 5: 158455 to 158486GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #5
Score
Expect
Identities
Gaps
Strand
54.5 bits(59)
0.015
31/32(97%)
0/32(0%)
Plus/Plus
Query  604     TAACTAAACTAATGCATCTGCTCAAAGGCATC  635
               |||||||||||||||||||||||| |||||||
Sbjct  158455  TAACTAAACTAATGCATCTGCTCAGAGGCATC  158486


Range 6: 284599 to 284630GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #6
Score
Expect
Identities
Gaps
Strand
54.5 bits(59)
0.015
31/32(97%)
0/32(0%)
Plus/Plus
Query  604     TAACTAAACTAATGCATCTGCTCAAAGGCATC  635
               |||||||||||||||||||||||| |||||||
Sbjct  284599  TAACTAAACTAATGCATCTGCTCAGAGGCATC  284630


Range 7: 218237 to 218305GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #7
Score
Expect
Identities
Gaps
Strand
48.2 bits(52)
0.62
54/70(77%)
2/70(2%)
Plus/Plus
Query  535     GTGGTGCGGAAGTTGATTAATCCAACCAAACGTGCC-TGCTTGTCTCTAGATAATACAAA  593
               ||||| |  |||||||||||||||| | ||| |||| ||| | |||  ||  ||| ||||
Sbjct  218237  GTGGTTCCCAAGTTGATTAATCCAAGC-AACTTGCCATGCCTATCTGGAGTCAATGCAAA  218295

Query  594     CTAAGATGTC  603
               || || ||||
Sbjct  218296  CTGAGTTGTC  218305


Range 8: 31205 to 31231GenBankGraphics
Next Match
Previous Match
First Match
Alignment statistics for match #8
Score
Expect
Identities
Gaps
Strand
45.5 bits(49)
7.6
26/27(96%)
0/27(0%)
Plus/Plus
Query  609    AAACTAATGCATCTGCTCAAAGGCATC  635
              ||||||||||||||||||| |||||||
Sbjct  31205  AAACTAATGCATCTGCTCAGAGGCATC  31231


Download
GenBank
Graphics
Sort by:  
Next
Previous
Descriptions
Triticum urartu cultivar G1812 clone BAC 288D18 chromosome 3AL, complete sequence 
Sequence ID: KC816724.1Length: 107101Number of Matches: 4
Range 1: 544 to 828GenBankGraphics
Next Match
Previous Match
Alignment statistics for match #1
Score
Expect
Identities
Gaps
Strand
142 bits(157)
4e-29
204/286(71%)
20/286(6%)
Plus/Plus
Query  433  GAAGTTCGCAGCCACATTGTCGTTGTGGAGTCTATATCAACCACCATTCTTAAAATTAGC  492
            |||||| || ||   | | | | || |||||||| |||||||||| ||||||||||||||
Sbjct  544  GAAGTTAGCTGCGGGAATATAGCTGCGGAGTCTAGATCAACCACCGTTCTTAAAATTAGC  603

Query  493  TACGTGTGATCAACTTGGTCCACCACCTACGCAAGGAGGTGCGTGGTGCGGAAGTTGATT  552
             |  |||| || | ||| ||||| ||  ||||||||||||||||||| |  |||||||||
Sbjct  604  AATATGTGGTC-ATTTGCTCCACTACGGACGCAAGGAGGTGCGTGGTTCCAAAGTTGATT  662

Query  553  -AATCCAAC--CAAACGTGCCTGCTTGTCTCTAGATAATAC-AAACTAAGATGTCTAACT  608
             ||||||||  | || |   |||   || |||||| ||| | ||||| || | | |||||
Sbjct  663  AAATCCAACGACCAAGGCAGCTGTGGGTTTCTAGACAATGCAAAACTTAGCTATTTAACT  722

Query  609  AAACTAATGCATCTGCTCA---------------AAGGCATCAATCGCATCGGCGCAAGC  653
            | ||||||||| ||| | |                || |||| || |||||  | |||| 
Sbjct  723  ACACTAATGCACCTGTTTAGAGAGTTGGGGCATCGAGCCATCGATGGCATCATCTCAAGT  782

Query  654  GTTCAGAGGTCAGACGACCATTGGCAAGCACCAAAGGACGGGGCCA  699
            || | || | |    | | |||||||||||||||||||| ||||||
Sbjct  783  GTCCTGACGCCCATTGTCTATTGGCAAGCACCAAAGGACAGGGCCA  828"""

    # Intentionally a nonstandard format.
    atgsp_example_results = """Query= 
Length=816
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

lcl|Chr3                                                              1480    0.0  
lcl|Chr1                                                               167    1e-38
lcl|Chr5                                                               163    2e-37
lcl|Chr2                                                               161    7e-37
lcl|Super_scaffold_579                                                 159    2e-36
lcl|Chr6                                                               145    7e-32
lcl|Chr7                                                               135    4e-29
lcl|Chr4                                                               110    2e-21


>Chr3    Length=627182665

 Score = 1480 bits (801),  Expect = 0.0
 Identities = 812/817 (99%), Gaps = 2/817 (0%)
 Strand=Plus/Minus

Query  1          TGCGATGACGGAAAAAAAAA-GGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTA  59
                  ||||||||| |||||||||| |||| ||||||||||||||||||||||||||||||||||
Sbjct  585890008  TGCGATGAC-GAAAAAAAAACGGTGATGGGAGTATGACGAAAATAAACCAGCGAAAATTA  585889950

Query  60         ACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTG  119
                  ||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889949  ACCGCGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGATTTTCACAACCCAATTTG  585889890

Query  120        CCTATTTTATTTATTCATTGGAAAACTTGAAGAATGGAAATTTTGTGGTGCGATAATAGC  179
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889889  CCTATTTTATTTATTCATTGGAAAACTTGAAGAATGGAAATTTTGTGGTGCGATAATAGC  585889830

Query  180        ATATTGTGTTGATGTTGTCAAATTTCTTCCCACCAACCGTGCCCTCCTCCAGTTGATTGT  239
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889829  ATATTGTGTTGATGTTGTCAAATTTCTTCCCACCAACCGTGCCCTCCTCCAGTTGATTGT  585889770

Query  240        ATGCACAAGCGTGCCACCCCATGACGGTAGGAAACGAGCACAACTGCACCATTGTTTAAA  299
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889769  ATGCACAAGCGTGCCACCCCATGACGGTAGGAAACGAGCACAACTGCACCATTGTTTAAA  585889710

Query  300        TCATCTCGTCAGAGCTTCATCAGTTTATAATGGCATATGTTTGTCCACAAGCGCATCACC  359
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889709  TCATCTCGTCAGAGCTTCATCAGTTTATAATGGCATATGTTTGTCCACAAGCGCATCACC  585889650

Query  360        GGCCCCGGTGCGCCAAATATATTTTGTAACCACCCTTATCCTTTCAATCAACCATCACGA  419
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889649  GGCCCCGGTGCGCCAAATATATTTTGTAACCACCCTTATCCTTTCAATCAACCATCACGA  585889590

Query  420        GGAAATAGTCTGAGAAGTTCGCAGCCACATTGTCGTTGTGGAGTCTATATCAACCACCAT  479
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889589  GGAAATAGTCTGAGAAGTTCGCAGCCACATTGTCGTTGTGGAGTCTATATCAACCACCAT  585889530

Query  480        TCTTAAAATTAGCTACGTGTGATCAACTTGGTCCACCACCTACGCAAGGAGGTGCGTGGT  539
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889529  TCTTAAAATTAGCTACGTGTGATCAACTTGGTCCACCACCTACGCAAGGAGGTGCGTGGT  585889470

Query  540        GCGGAAGTTGATTAATCCAACCAAACGTGCCTGCTTGTCTCTAGATAATACAAACTAAGA  599
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889469  GCGGAAGTTGATTAATCCAACCAAACGTGCCTGCTTGTCTCTAGATAATACAAACTAAGA  585889410

Query  600        TGTCTAACTAAACTAATGCATCTGCTCAAAGGCATCAATCGCATCGGCGCAAGCGTTCAG  659
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889409  TGTCTAACTAAACTAATGCATCTGCTCAAAGGCATCAATCGCATCGGCGCAAGCGTTCAG  585889350

Query  660        AGGTCAGACGACCATTGGCAAGCACCAAAGGACGGGGCCATACGCTTAGAGACACCTACA  719
                  |||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||
Sbjct  585889349  AGGTCAGACGACCATTGGCAAGCACCAAAGGACGGGGCCATACGCTTAGAAACACCTACA  585889290

Query  720        TGGCCGGGTATTAATTACTAGAGTGTTGCAACGCGTTCGCAAAGTCAGCACACGGAGATT  779
                  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  585889289  TGGCCGGGTATTAATTACTAGAGTGTTGCAACGCGTTCGCAAAGTCAGCACACGGAGATT  585889230

Query  780        GGTGCCGCAATGTAATTTGACTCTGGCATGACAGTCC  816
                  |||||||||||||||||||||||||||||||||||||
Sbjct  585889229  GGTGCCGCAATGTAATTTGACTCTGGCATGACAGTCC  585889193


 Score =  169 bits (91),  Expect = 4e-39
 Identities = 100/104 (96%), Gaps = 2/104 (2%)
 Strand=Plus/Minus

Query  1          TGCGATGACGGAAAAAAAAA-GGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTA  59
                  ||||||||| |||||||||| || ||||||||||||||||||||||||||||||||||||
Sbjct  148076245  TGCGATGAC-GAAAAAAAAACGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTA  148076187

Query  60         ACCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGATT  103
                  ||||||||||||||||| ||||||||||||||||||||||||||
Sbjct  148076186  ACCGCGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGATT  148076143


 Score =  161 bits (87),  Expect = 7e-37
 Identities = 96/100 (96%), Gaps = 1/100 (1%)
 Strand=Plus/Plus

Query  3          CGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACC  62
                  ||||||||| ||||||| || |||||||||||||||||||||||||||||||||||||||
Sbjct  396782493  CGATGACGG-AAAAAAACGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACC  396782551

Query  63         GCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGAT  102
                  |||||||||||||| |||||||||||||||||||||||||
Sbjct  396782552  GCGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  396782591


 Score =  159 bits (86),  Expect = 2e-36
 Identities = 97/102 (95%), Gaps = 1/102 (1%)
 Strand=Plus/Minus

Query  1          TGCGATGACGGAAAAAAAAAGGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  60
                  ||||||||  |||||||||||| |||||||||||||||||||||||||||||||||||||
Sbjct  528388742  TGCGATGA-TGAAAAAAAAAGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAA  528388684

Query  61         CCGCGCGGACTATTCATCAAGTCGTCCATTAGGAGTAGAGAT  102
                  ||| |||||||||||| |||||||||||||||||||||||||
Sbjct  528388683  CCGTGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT  528388642


 Score =  156 bits (84),  Expect = 3e-35
 Identit"""


    def test_ncbi_length(self):
        parser = PairwiseParser()
        hsps = parser.parse(self.ncbi_test_results)
        self.assertEqual(len(hsps), 9)

    def test_ncbi_query_seq_6(self):
        parser = PairwiseParser()
        hsps = parser.parse(self.ncbi_test_results)
        self.assertEqual(hsps[6].q_seq, "GTGGTGCGGAAGTTGATTAATCCAACCAAACGTGCC-TGCTTGTCTCTAGATAATACAAACTAAGATGTC")

    def test_ncbi_midline_1(self):
        parser = PairwiseParser()
        hsps = parser.parse(self.ncbi_test_results)
        self.assertEqual(hsps[1].midline, "|||||||||||||| ||  |||| ||||||    || | |||| | ||||||| ||||| ||||||| |||||||||||||||||||||| |||       |||  | || | |     ||||||  || | ||      | |||||   ||| ||  |   | ||||| ||   ||  ||||||| |||  | || | ||  |||| ||||||||| ||||||| ||||||||||||||||||| |   || | | |||  ||||||||||||||| ||| |||||")

    def test_atgsp_length(self):
        parser = PairwiseParser()
        hsps = parser.parse(self.atgsp_example_results)
        self.assertEqual(len(hsps), 4)

    def test_atgsp_subject_seq_3(self):
        parser = PairwiseParser()
        hsps = parser.parse(self.atgsp_example_results)
        self.assertEqual(hsps[3].s_seq, "TGCGATGA-TGAAAAAAAAAGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTAACCGTGCGGACTATTCACCAAGTCGTCCATTAGGAGTAGAGAT")

if __name__ == "__main__":
    unittest.main()
