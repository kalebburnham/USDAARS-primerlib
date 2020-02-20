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

    Score =  569 bits (308),  Expect = 4e-159
    Identities = 478/561 (85%), Gaps = 8/561 (1%)
    Strand=Plus/Plus

    Query  2643       TATCTTCATTGTATTGATTTTATAACCGATTCCAAAATGTATTCTTAAAGGTACATCATC  2702
                      |||| ||||||||||||| |  ||||  |||   ||||| |||| ||||||||| |||||
    Sbjct  586624978  TATCCTCATTGTATTGATCTATTAACTAATTATTAAATGCATTCATAAAGGTACCTCATC  586625037

    Query  2703       GTAATTGATGATATATGGGATGAAAAAGTGTGGGAATTCATTAA-T-TGCGCTTTCTCCA  2760
                      ||||| |||||||||||| |||||||||  ||||| ||  |||| | ||| |||| ||||
    Sbjct  586625038  GTAATCGATGATATATGGAATGAAAAAGCATGGGAGTTACTTAAGTGTGC-CTTT-TCCA  586625095

    Score =  169 bits (91),  Expect = 2e-38
    Identities = 100/104 (96%), Gaps = 2/104 (2%)
    Strand=Plus/Minus


    Query  1          TGCGATGACGGAAAAAAAAA-GGTGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTA  59
                      ||||||||| |||||||||| || ||||||||||||||||||||||||||||||||||||
    Sbjct  148076245  TGCGATGAC-GAAAAAAAAACGGCGGTGGGAGTATGACGAAAATAAACCAGCGAAAATTA  148076187

    Note that it is not recommended to use this format.
    """

    def __init__(self):
        pass

    def parse(self, data):
        """
        Returns a list of the subject sequence, aka hit sequences (hseqs).
        """
        data = self.tokenize(data)

        hseqs = []
        hit_sequence = ''

        idx = 0
        while idx < len(data):
            if data[idx] == 'Score':
                # A new alignment is starting. Save the current
                # hit sequence if it is nonempty.
                if hit_sequence:
                    hseqs.append(hit_sequence)
                    hit_sequence = ''

            if data[idx] == 'Sbjct':
                # Sbjct line found. The sequence is 2 tokens away.
                idx += 2
                hit_sequence += data[idx].replace('-', '')

            idx += 1

        if hit_sequence:
            hseqs.append(hit_sequence)

        return hseqs

    def tokenize(self, data):
        """ Splits the data at every new line and whitespace and returns
        the list of tokens. """
        return data.split()
