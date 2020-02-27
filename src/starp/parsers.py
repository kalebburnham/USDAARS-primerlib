"""
License information goes here.
"""

import re
import regex
from bs4 import BeautifulSoup

from .exceptions import StarpError
from .models import Sequence, Snp

def get_parser(data):
    """
    Given the user's input data, return the correct Parser.

    It is highly recommended (by NCBI) not to use BLAST data pairwise
    data in other programs. That is only meant to be human-readable and
    is subject to significant changes in the future.

    All accepted formats MUST include a snps() function that returns the
    SNPs in HGVS's standard representation.

    Args:
        data: The user's input data.

    Returns:
        The correct parser for their data format.

    Raises:
        StarpError: if format is not recognized.
    """

    # If there are only two lines and the first begins with a '>', the
    # data is in FASTA format. Remove the first line to get the
    # sequence.
    if len(data.splitlines()) == 2:
        if data.startswith('>'):
            data = data.splitlines()[1]

    # Test for SnpSequence
    pattern = regex.compile(r'\w|\[.\/.\]')
    matched_chars = ''.join(regex.findall(pattern, data))
    if matched_chars == data:
        return SnpSequence(data)

    # Test for TwoAlleles
    lines = data.splitlines()
    if len(lines) == 4 and lines[0].startswith('>') and lines[2].startswith('>'):
        return TwoAlleles(data)

    # Test for Single Blast Sequence
    if '|' in data:
        return SingleBlastParser(data)

    # Format not recognized.
    raise StarpError("SNP Format Not Recognized")

class SnpSequence:
    """
    Defined as a single sequence of nucleotides with SNPs specified
    by [a/b], a != b.

    Accepted alphabet is {A, C, G, T, -, /, [, ]}

    The accepted Grammar of this format is:
    S -> nS | [a/b]S | empty_string
    where n, a, b in {A, C, G, T, -} and a != b.

    """

    def __init__(self, data):
        self.data = data

        ALPHABET = 'ACGT-/[]'
        SNP_ALPHABET = 'ACGT-'
        SNP_PATTERN = r'\[[' + SNP_ALPHABET + r']\/[' + SNP_ALPHABET + r']\]'

        """ Search for invalid characters. """
        if regex.search(r'[^ACGT\-/\[\]]', data):
            raise StarpError(('The accepted alphabet for a SNP sequence is '
                              '{A, C, G, T, -, /, [, ]}.'))

        """ Search for SNPs with the same characters, eg [A/A] and [-/-]. """
        if regex.search(r'\[([ACGT-])\/(\1)\]', data):
            raise StarpError('SNPs must contain different letters from '
                             '' + SNP_ALPHABET + '!')

        """Search for anything that doesn't match the grammar.
        Kind of a hack, but this is done by searching for the
        grammar that's wanted. Then if there are unmatched chars,
        the length of matched_chars is less than the length of the
        data. """
        pattern = regex.compile('[' + SNP_ALPHABET + ']|' + SNP_PATTERN)
        matched_chars = ''.join(regex.findall(pattern, data))
        if len(matched_chars) < len(data):
            raise StarpError("Invalid syntax.")

        self.allele1_aligned = ''
        self.allele2_aligned = ''
        self.allele1 = ''
        self.allele2 = ''

        # Call snps() to set the four preceding attributes.
        self.snps()

    def tokenize(self, data):
        # Split string into list of chars.
        return list(data)

    def snps(self):
        """
        Parses and returns the SNP objects from the data.
        """
        data = self.tokenize(self.data)

        self.allele1_aligned = []  # Will be joined into string later.
        self.allele2_aligned = []
        snps = []

        idx = 0 # Current index on the SNP sequence.
        pos = 0 # Position relative to the first allele
        while idx < len(data):
            c = data[idx]
            if c in 'ACGTN':
                self.allele1_aligned.append(c)
                self.allele2_aligned.append(c)
            elif c == '[':
                self.allele1_aligned.append(data[idx+1])
                self.allele2_aligned.append(data[idx+3])
                if data[idx+1] == '-':
                    # Insertion SNP
                    descriptor = f'.{pos}ins{data[idx+3]}'
                    snp = Snp(descriptor)
                elif data[idx+3] == '-':
                    # Deletion SNP
                    descriptor = f'.{pos}del'
                    snp = Snp(descriptor)
                    snp.ref_nucleotide = data[idx+1]
                else:
                    # Substitution SNP
                    descriptor = f'.{pos}{data[idx+1]}>{data[idx+3]}'
                    snp = Snp(descriptor)
                pos += 1

                snps.append(snp)

                # Place the idx on the ']' so the next increment reads the next token.
                idx += 4
            else:
                raise StarpError(('Invalid characters in Snp Sequence. The accepted alphabet '
                                  'is {A, C, G, T, -, /, [, ]}.'))

            idx += 1

        self.allele1_aligned = ''.join(self.allele1_aligned)
        self.allele2_aligned = ''.join(self.allele2_aligned)
        self.allele1 = self.allele1_aligned.replace('-', '')
        self.allele2 = self.allele2_aligned.replace('-', '')

        return snps

class TwoAlleles:
    """
    Parses SNPs out of two alleles entered as

    >Allele 1
    NNNNNNNNNNANNNNNNNNNN
    >Allele 2
    NNNNNNNNNNBNNNNNNNNNN

    Each allele must be a continuous string (no new lines).

    The only accepted alphabet for the alleles {A, C, G, T, -, N}

    Attributes:
        allele1_aligned: The original Allele 1 data with '-' as
            placeholders for insertions and deletions.
        allele2_aligned: The original Allele 2 data with '-' as
            placeholders for insertions and deletions.
        allele1: The first allele but without '-'s.
        allele2: The second allele but without '-'s.
    """
    def __init__(self, data):
        """
        Args:
            data: A string of four lines in the following format.

                    >Allele 1
                    NNNNNNNNNNANNNNNNNNNN
                    >Allele 2
                    NNNNNNNNNNBNNNNNNNNNN

        Raises:
            StarpError: The input alleles are of differing lengths.
            StarpError: The alleles have characters not in "ACGTN-"
        """
        lines = data.splitlines()
        self.allele1_aligned = Sequence(lines[1].rstrip())
        self.allele2_aligned = Sequence(lines[3].rstrip())

        self.allele1 = Sequence(str(self.allele1_aligned).replace('-', ''))
        self.allele2 = Sequence(str(self.allele2_aligned).replace('-', ''))

        if len(self.allele1_aligned) != len(self.allele2_aligned):
            raise StarpError("The alleles are of differing lengths.")

        if regex.search('[^ACGTN-]', str(self.allele1_aligned)):
            raise StarpError("Allele 1 has invalid characters.\
                    The only valid characters are A, C, G, T, N, and -")

        if regex.search('[^ACGTN-]', str(self.allele2_aligned)):
            raise StarpError("Allele 2 has invalid characters.\
                    The only valid characters are A, C, G, T, N, and -")

    def snps(self):
        """
        Computes the SNPs between the two aligned alleles.

        Args:
            None

        Returns:
            A list of Snp objects representing all Snps between the
            aligned allele sequences.

        Raises:
            None
        """
        snps = list()
        position = 0

        # Insertions throw off the positioning of Snps, so we need to
        # keep track of them and offset the positions when there are
        # many insertions.

        # Example: Allele1 = TGACACGTACGT
        # Insertion of 'A' at position 2 and 'G' at position 8 gives
        # TAGACACGGTACGT.
        num_insertions = 0

        # Constructs SNP objects by comparing element-wise the characters in
        # the aligned sequences.
        for pair in zip(str(self.allele1_aligned), str(self.allele2_aligned)):
            if pair[0] == pair[1]:
                pass # No snp here.
            elif pair[0] == 'N' or pair[1] == 'N':
                # No SNP will have an 'N'.
                pass
            elif pair[0] == '-':
                # Insertion SNP
                descriptor = '.' + str(position) + 'ins' + pair[1]
                snps.append(Snp(descriptor))
                num_insertions += 1
            elif pair[1] == '-':
                # Deletion SNP. The descriptor does not contain any
                # information on reference/new nucleotides, but the
                # reference nucleotide is known so it is manually
                # entered.
                descriptor = '.' + str(position) + 'del'
                snp = Snp(descriptor)
                snp.ref_nucleotide = str(self.allele1_aligned[position])
                snps.append(snp)
            else:
                # Substitution SNP
                descriptor = '.' + str(position) + str(pair[0]) + '>' + str(pair[1])
                snps.append(Snp(descriptor))
            position += 1

        return snps

class SingleBlastParser:
    """ Parses a single BLAST alignment.
    For example,

    Query  2643       TATCTTCATTGTATTGATTTTATAACCGATTCCAAAATGTATTCTTAAAGGTACATCATC  2702
                      |||| ||||||||||||| |  ||||  |||   ||||| |||| ||||||||| |||||
    Sbjct  586624978  TATCCTCATTGTATTGATCTATTAACTAATTATTAAATGCATTCATAAAGGTACCTCATC  586625037

    Query  2703       GTAATTGATGATATATGGGATGAAAAAGTGTGGGAATTCATTAA-T-TGCGCTTTCTCCA  2760
                      ||||| |||||||||||| |||||||||  ||||| ||  |||| | ||| |||| ||||
    Sbjct  586625038  GTAATCGATGATATATGGAATGAAAAAGCATGGGAGTTACTTAAGTGTGC-CTTT-TCCA  586625095

    is valid but

    Score =  569 bits (308),  Expect = 4e-159
    Identities = 478/561 (85%), Gaps = 8/561 (1%)
    Strand=Plus/Plus

    Query  2643       TATCTTCATTGTATTGATTTTATAACCGATTCCAAAATGTATTCTTAAAGGTACATCATC  2702
                      |||| ||||||||||||| |  ||||  |||   ||||| |||| ||||||||| |||||
    Sbjct  586624978  TATCCTCATTGTATTGATCTATTAACTAATTATTAAATGCATTCATAAAGGTACCTCATC  586625037

    Query  2703       GTAATTGATGATATATGGGATGAAAAAGTGTGGGAATTCATTAA-T-TGCGCTTTCTCCA  2760
                      ||||| |||||||||||| |||||||||  ||||| ||  |||| | ||| |||| ||||
    Sbjct  586625038  GTAATCGATGATATATGGAATGAAAAAGCATGGGAGTTACTTAAGTGTGC-CTTT-TCCA  586625095

    is not.
    """

    def __init__(self, data):
        self.data = data
        self.allele1 = ''
        self.allele2 = ''
        self.allele1_aligned = ''
        self.allele2_aligned = ''

        # Call snps() to set the 4 latter attributes.
        self.snps()

    def snps(self):
        tokens = self.tokenize(self.data)
        query = ''
        sbjct = ''

        i = 0
        while i < len(tokens):
            if tokens[i].upper() == 'QUERY':
                # Query subsequence is 2 tokens away
                i += 2
                query += tokens[i]
            elif tokens[i].upper() == 'SBJCT':
                # Sbjct subsequence is 2 tokens away
                i += 2
                sbjct += tokens[i]
            i += 1
        
        # Put the query and sbjct sequences in the proper format and
        # utilize the Two Alleles Parser.
        two_alleles = TwoAlleles(f'>Query\n{query}\n>Sbjct\n{sbjct}')
        self.allele1 = two_alleles.allele1
        self.allele2 = two_alleles.allele2
        self.allele1_aligned = two_alleles.allele1_aligned
        self.allele2_aligned = two_alleles.allele2_aligned

        return two_alleles.snps()

    def tokenize(self, data):
        """ Splits the data at every new line and whitespace and returns
        the list of tokens. """
        tokens = data.split()
        for token in tokens:
            token.replace('|', '')
        return tokens

# *********
# Nontarget parsers
class XmlParser:
    """ Returns a list of hit sequences from BLAST XML data. """

    def __init__(self):
        pass

    def parse(self, xml):
        """ Returns a list of hit sequences from the xml data. """
        soup = BeautifulSoup(xml, 'xml')
        hseqs = soup.find_all('Hsp_hseq')
        return [hseq.get_text().replace('-', '') for hseq in hseqs]

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
