"""
License information goes here.
"""

import re
import textwrap

from bs4 import BeautifulSoup

from .models import Hsp

# Any parsers needed for parsing BLAST data should be placed here.
# Note that it is not recommended to parse the text output from BLAST
# since that format has changed in the past and could change in the
# future, making our program unusable until the parsers are updated.

class XmlParser:
    """ Take XML data from BLAST searches and put them into HSP
    objects. Works with BLAST XML version 1.0.
    """

    def __init__(self, size=80):
        """ size is the length of each line in a human-readable
        pairwise output. Not actually used in the parsing.
        """
        self.size = size

    def parse(self, xml):
        """ Returns a list of hit sequences from the xml data. """
        soup = BeautifulSoup(xml, 'xml')
        hseqs = soup.find_all('Hsp_hseq')
        return [hseq.get_text() for hseq in hseqs]

class PairwiseParser:
    """
    Pairwise data is the human-readable blast output. It looks like

            100  TGACTGAGTGTGTACCAGATAGCTTGCGTCAGTCGAGTCAT 141
                 |||||||||||||||||||||||||||||||||||||||||
            4930 TGACTGAGTGTGTACCAGATAGCTTGCGTCAGTCGAGTCAT 4971

            353   TGTCTGTACAT....
                  |||||||||||....
            23490 TGTCTGTACAT....

    Note that it is not recommended to use this format.
    """

    def __init__(self):
        pass

    def parse(self, data):
        """
        data - the pairwise output from BLAST.
        Since each site uses different formats, we have to be careful
        and only look for the query, midline, and subjects from each
        result and must ignore all the metadata of each hit. This is
        achieved by reading line-by-line.

        Keep in mind it is NOT recommended to use this type of data
        in the primer programs because the outputs can change over time,
        they differ between BLAST implementations, and the user input
        is non-standardized. It exists solely because some sites don't
        offer other formats (See: Aegilops Tauschii Genome Sequencing
        Project, Oct. 2019).
        """
        hsps = list()

        lines = data.splitlines()
        iterator = iter(lines)

        query_lines = list()
        mid_lines = list()
        subject_lines = list()

        try:
            while True:
                line = next(iterator)
                if line.startswith('Query '):
                    preface_length = self.get_preface_length(line)
                    start, stop = self.get_seq_start_stop(line)
                    query_line = line
                    mid_line = next(iterator)
                    subject_line = next(iterator)
                    blank_line = next(iterator)

                    query_lines.append(query_line[start:stop])
                    mid_lines.append(mid_line[start:stop])
                    subject_lines.append(subject_line[start:stop])
                elif query_lines:
                    # Save the hit, if there is one.
                    query_seq = ''.join(query_lines)
                    midline = ''.join(mid_lines)
                    subject_seq = ''.join(subject_lines)
                    # Unfortunately, we have to ignore all the metadata
                    # because this format can be wildly different
                    # between sites. Also, how the user inputs this data
                    # is not standardized.
                    hsps.append(Hsp(hit_accession='',
                                    q_from=0,
                                    q_to=0,
                                    s_from=0,
                                    s_to=0,
                                    q_seq=query_seq,
                                    s_seq=subject_seq,
                                    midline=midline,
                                    web_alignment=''))

                    # Reset our arrays for the next hit.
                    query_lines = list()
                    mid_lines = list()
                    subject_lines = list()
        except StopIteration:
            # Save the hit, if there is one.
            if query_lines:
                query_seq = ''.join(query_lines)
                midline = ''.join(mid_lines)
                subject_seq = ''.join(subject_lines)
                hsps.append(Hsp(hit_accession='',
                                q_from=0,
                                q_to=0,
                                s_from=0,
                                s_to=0,
                                q_seq=query_seq,
                                s_seq=subject_seq,
                                midline=midline,
                                web_alignment=''))

        return hsps

    def get_preface_length(self, query_line):
        """
        Return the length of the preface of each pairwise line using
        a query line. Example: "Query  433    GAAGTTCGC...." contains
        5 characters for query, 3 characters for '433', and 6 spaces
        before we get to the actual sequence. So, preface length is 11.
        """
        pattern = re.compile('[^ACGT-]+', flags=re.IGNORECASE)
        match = re.match(pattern, query_line)
        return len(match.group(0))

    def get_seq_start_stop(self, query_line):
        """
        Return the start and stop of the sequence on each line, using the
        query line for reference.
        Example: "Query 1  GTGTAATAGC  10" would have start 9 and stop 19.
        """
        pattern = re.compile('[ACGT-]+', flags=re.IGNORECASE)
        match = re.search(pattern, query_line)
        return match.start(), match.end()

if __name__ == "__main__":
    data = """Triticum monococcum subsp. monococcum cultivar DV92 Sr35 region, genomic sequence
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
Sbjct  783  GTCCTGACGCCCATTGTCTATTGGCAAGCACCAAAGGACAGGGCCA  828
"""
    parser = PairwiseParser()
    hsps = parser.parse(data)
    print(hsps[8].q_seq)
    print(hsps[8].midline)
    print(hsps[8].s_seq)
    print(len(hsps))
