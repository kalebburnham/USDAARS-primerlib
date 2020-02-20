"""
License information goes here.
"""

import re
import regex
from bs4 import BeautifulSoup

from .exceptions import StarpError
from .models import Hsp, Sequence, Snp

ALPHABET = 'ACGT-/()'
SNP_ALPHABET = 'ACGT-'

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

    # Format not recognized.
    raise StarpError("SNP Format Not Recognized")

class SnpSequence:
    """
    Defined as a single sequence of nucleotides with SNPs specified
    by (a/b), a != b.

    Accepted alphabet is {A, C, G, T, -, /, (, )}

    The accepted Grammar of this format is:
    S -> nS | (a/b)S | empty_string
    where n, a, b in {A, C, G, T, -} and a != b.

    Raises ValueError if the data contains invalid characters.

    Raises SyntaxError if Snps contain the same character or the data
    has incorrect grammar.
    """
    def __init__(self, data):
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

        # Convert data into TwoAlleles format
        data_iter = iter(data)
        allele1 = list()
        allele2 = list()

        while True:
            try:
                c = next(data_iter)
            except StopIteration:
                break

            if c == '[':
                allele1.append(next(data_iter))
            elif c == '/':
                allele2.append(next(data_iter))
            elif c == ']':
                # Ignore this.
                pass
            else:
                # This character is common to both alleles.
                allele1.append(c)
                allele2.append(c)

        self.allele1_aligned = Sequence(''.join(allele1))
        self.allele2_aligned = Sequence(''.join(allele2))

        self.allele1 = Sequence(''.join(allele1).replace('-', ''))
        self.allele2 = Sequence(''.join(allele2).replace('-', ''))

        self._twoAllelesFormat = TwoAlleles('>Allele1\n'
                                            + str(self.allele1_aligned)
                                            + '\n>Allele2\n'
                                            + str(self.allele2_aligned))

    def snps(self):
        """
        Returns the SNP objects from the data.
        """
        return self._twoAllelesFormat.snps()

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
                descriptor = '.' + str(position-num_insertions) + 'ins' + pair[1]
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

# *********
# Nontarget parsers
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

    def parse2(self, xml_data):
        """ Return a list of strings of HSPs with web alignment that
        looks like this:

        GTAGTCATGCGTATGCGTATGCACATTAAGTCGTGTACTGACTG
        ||||||||||||||||||||||||  |||||||||||||| |||
        GTAGTCATGCGTATGCGTATGCACCGTAAGTCGTGTACTGACTG

        GTAGTCATGCGTATGCGTATGCACATTAAGTCGTGTACTGACTG
        ||||||||||||||||||||||||  |||||||||||||| |||
        GTAGTCATGCGTATGCGTATGCACCGTAAGTCGTGTACTGACTG

        Blast result should be in XML.
        """
        hsps = list()

        # The Internet Explorer viewer displays XML with
        # leading dashes to collapse sections.
        # If these dashes are copied, it causes issues with
        # the parsing. They are removed here.
        xml_lines = xml_data.splitlines()
        for line in xml_lines:
            if line.startswith('-'):
                line = line[1:]
        xml_data = ''.join(xml_lines)

        # Change the dashes in these examples to underscores,
        # since dashes can mess up the parsing.
        xml_data = re.sub(r'query-', 'query_', xml_data)
        xml_data = re.sub(r'hit-', 'hit_', xml_data)

        soup = BeautifulSoup(xml_data, 'lxml-xml')
        for hit in soup.find_all('Hit'):

            hit_accession = hit.Hit_accession.text
            for hsp in hit.find_all('Hsp'):
                if int(hsp.Hsp_hit_frame.text) == -1:
                    # If the hit frame is -1, the query string is
                    # reversed in the output
                    q_from = int(hsp.Hsp_query_to.text)
                    q_to = int(hsp.Hsp_query_from.text)
                    s_from = int(hsp.Hsp_hit_to.text)
                    s_to = int(hsp.Hsp_hit_from.text)

                    q_seq = hsp.Hsp_qseq.text[::-1]
                    s_seq = hsp.Hsp_hseq.text[::-1]
                    midline = hsp.Hsp_midline.text[::-1]
                else:
                    q_from = int(hsp.Hsp_query_from.text)
                    q_to = int(hsp.Hsp_query_to.text)
                    s_from = int(hsp.Hsp_hit_from.text)
                    s_to = int(hsp.Hsp_hit_to.text)

                    q_seq = hsp.Hsp_qseq.text
                    s_seq = hsp.Hsp_hseq.text
                    midline = hsp.Hsp_midline.text

                q = self._str2list(q_seq)
                s = self._str2list(s_seq)
                m = self._str2list(midline)

                rows = list()
                for i in range(len(q)):
                    rows.append(q[i] + '\n' + m[i] + '\n' + s[i] + '\n\n')

                web_alignment = ''.join(rows)

                hsps.append(Hsp(hit_accession, q_from, q_to, s_from, s_to,
                                q_seq, s_seq, midline, web_alignment))
        if hsps:
            if hsps[0].midline.count('|') == len(hsps[0].q_seq):
                del hsps[0]

        return hsps

    def _str2list(self, s):
        """ 
        Break s into a list of lines where each line is at most
        self.size characters long.
        """
        lines = list()
        full_rows = int(len(s) / self.size)
        for i in range(full_rows):
            lines.append(s[i*self.size:(i+1)*self.size])
        lines.append(s[full_rows*self.size:])
        return lines

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