import math
import re
import json

import nestedloop.utils as utils

class Additive:
    """
    A class that automatically calculates PCR conditions given a
    target sequence upon initialization.

    Refer to docs/NL-PCR_reaction_conditions.docx

    Attributes:
        sequence: The user-defined target sequence.
        max_gc: The maximum gc percentage of any 100bp segment
        min_gc: The minimum gc percentage of any 100bp segment.
        additive: The amount of additive needed in PCR.
        pcr_temperatures: The temperature sequence for PCR.
        n: Unknown. Refer to docs/NL-PCR_reaction_conditions.docx
    """

    def __init__(self, target_sequence):
        self.sequence = target_sequence
        self.max_gc = self._max_gc()
        self.min_gc = self._min_gc()
        self.additive = self._additive()
        self.pcr_temperatures = self._pcr_temperatures()

        # What is n? It is on the PCR conditions document.
        self.n = 1 + math.ceil((len(target_sequence) / 3000) - 0.3)

    def _additive(self):
        """ Units are in ul. """
        additive = None

        if math.isnan(self.max_gc):
            return float('nan')

        if self.max_gc >= 0.95:
            additive = "6.5-9.0"
        elif self.max_gc >= 0.90:
            additive = "6.0-8.5"
        elif self.max_gc >= 0.85:
            additive = "5.5-8.0"
        elif self.max_gc >= 0.80:
            additive = "3.5-6.0"
        elif self.max_gc >= 0.75:
            additive = "3.0-5.5"
        elif self.max_gc >= 0:
            additive = "2.5-5.0"

        return additive

    def _max_gc(self):
        """ Returns the max GC content of any 100 bp subsequence
        in sequence. """

        # What to do about Ns?
        if len(self.sequence) < 100:
            # What to do?
            return float('nan')

        # Get all 100 bp sequence from the target sequence.
        fragments = [self.sequence[i:i+100] for i in range(len(self.sequence)-99)]
        gcs = [fragment.gc for fragment in fragments]
        return max(gcs)

    def _min_gc(self):
        if len(self.sequence) < 100:
            return float('nan')

        fragments = [self.sequence[i:i+100] for i in range(len(self.sequence)-99)]
        gcs = [fragment.gc for fragment in fragments]
        return min(gcs)

    def _pcr_temperatures(self):
        """ Returns temperatures in a pre-formatted HTML form."""
        temperatures = None

        if self.max_gc >= 0.95:
            if self.min_gc <= 0.10:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 58/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 59/30\"\n"
            else:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
        elif self.max_gc >= 0.90:
            if self.min_gc <= 0.10:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 59/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
            else:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 61/30\"\n"
        elif self.max_gc >= 0.85:
            if self.min_gc <= 0.10:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 61/30\"\n"
            else:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n"
        elif self.max_gc >= 0.80:
            temperatures = "1) 69/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n"
        else:
            temperatures = "1) 68/30\"\n2) 66/30\"\n3) 64/30\"\n4) 62/30\"\n"
            #temperatures = "68/30\"; 66/30\"; 64/30\"; 62/30\";"

        return temperatures

class Hsp:
    """ Represents a High Scoring Pair. These are output from BLAST in
    other formats. So, we need parsers to convert them to this format.
    More attributes can be needed, but these are what are needed so
    far. """

    def __init__(self, hit_accession, q_from, q_to, s_from, s_to, q_seq, s_seq,
                 midline, web_alignment):
        self.hit_accession = hit_accession
        self.q_from = q_from
        self.q_to = q_to
        self.s_from = s_from
        self.s_to = s_to
        self.q_seq = Sequence(q_seq)
        self.s_seq = Sequence(s_seq)
        self.midline = midline
        self.web_alignment = web_alignment

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and
            self.hit_accession == other.hit_accession and
            self.q_from == other.q_from and
            self.q_to == other.q_to and
            self.s_from == other.s_from and
            self.s_to == other.s_to and
            self.q_seq == other.q_seq and
            self.s_seq == other.s_seq and
            self.midline == other.midline and
            self.web_alignment == other.web_alignment
        )

    def __hash__(self):
        return hash((self.hit_accession, self.q_from, self.q_to, self.s_from,
                     self.s_to, self.q_seq, self.s_seq, self.midline,
                     self.web_alignment))

    def __str__(self):
        output = []
        output.append('<Hit-accession>' + self.hit_accession + '</Hit-accession>\n')
        output.append('    <q_from>' + str(self.q_from) + '</q_from>\n')
        output.append('    <q_to>' + str(self.q_to) + '</q_to>\n')
        output.append('    <s_from>' + str(self.s_from) + '</s_from>\n')
        output.append('    <s_to>' + str(self.s_to) + '</s_to>\n')
        output.append('    <midline>' + self.midline + '</midline>\n')
        output.append(re.sub('<br>', '\\n', '' + self.web_alignment))
        return ''.join(output)

class Pair:
    """ Container for pairs of forward and reverse primers. It is
    assumed the forward primer is on the plus strand and the reverse
    primer is on the minus strand, both written 5' -> 3' on the
    corresponding strand.

    Attributes:
        forward_primer: A Primer object on the plus strand.
        reverse_primer: A Primer object on the minus strand.
        complementary_score: The maximum number of complementary
            nucleotides between these two pairs.
        contig_complementary_score: The maximum number of contiguous
            complementary nucleotides between these two pairs.
        identity: The maximum number of identical nucleotides between
            these two pairs.
        contig_identity: The maximum number of contiguous identical
            nucleotides between these two pairs.
        distance: The length of the region spanned by these primers.

    Public Methods:
        additive: Returns an additive object containing the PCR
            condition information.
     """

    def __init__(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self._additive = None
        self._complementary_score = None
        self._contig_complementary_score = None
        self._identity = None
        self._contig_identity = None

    # We add the +1 to the starting position to make it one-indexed.
    def __str__(self):
        return ('Forward: ' + self.forward_primer.sequence
                + ' Reverse: ' + self.reverse_primer.sequence
                + ' Region: ' + str(self.forward_primer.start+1) + '-'
                + str(self.reverse_primer.end))

    def __repr__(self):
        return (f'Pair(forward={self.forward_primer.sequence}, '
                f'reverse={self.reverse_primer.sequence}, '
                f'span=({self.forward_primer.start}-{self.reverse_primer.end}), '
                f'complementary_score={self.complementary_score}, '
                f'contig_complementary_score={self.contig_complementary_score})')

    def toJSON(self):
        return {'forward_primer' : {'sequence': self.forward_primer.sequence,
                                               'span': self.forward_primer.span,
                                               'tm': self.forward_primer.tm,
                                               'gc': self.forward_primer.gc},
                           'reverse_primer': {'sequence': self.reverse_primer.sequence,
                                               'span': self.reverse_primer.span,
                                               'tm': self.reverse_primer.tm,
                                               'gc': self.reverse_primer.gc},
                           'distance': self.distance,
                           'additive': self.additive().additive,
                           'pcr_temperatures': self.additive().pcr_temperatures,
                           'n': self.additive().n}

    def additive(self, ref_sequence=None):
        """ Returns an Additive object with target sequence between
        the primers in this pair. Need the reference sequence to find
        the target sequence.
        """
        if self._additive is not None:
            return self._additive

        self._additive = Additive(ref_sequence[self.forward_primer.start:self.reverse_primer.end])
        return self._additive

    @property
    def complementary_score(self):
        """ See the Sequence implementation of this method.
        This assumes the reverse primer is on the -1 strand. """
        if self._complementary_score is not None:
            return self._complementary_score

        self._complementary_score = utils.complementary_score(
                                        self.forward_primer,
                                        self.reverse_primer.reverse())

        return self._complementary_score

    @property
    def contig_complementary_score(self):
        """ See the Sequence implementation of this method.
        This assumes the reverse primer is on the -1 strand. """
        if self._contig_complementary_score is not None:
            return self._contig_complementary_score

        self._contig_complementary_score = utils.contig_complementary_score(
                                               self.forward_primer,
                                               self.reverse_primer.reverse())

        return self._contig_complementary_score

    @property
    def distance(self):
        """ Returns the length of the amplicon formed by these
        primers. """
        return self.reverse_primer.end - self.forward_primer.start

    @property
    def identity(self):
        """ Calculate the identity between the pair sequences, where
        the identity is the number of pairs they have in common in a
        maximal alignment. Similar to complementary score, except we
        check for sameness in the original sequences."""
        if self._identity is not None:
            return self._identity

        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -1000
        self._identity = aligner.score(self.forward_primer.sequence,
                                       self.reverse_primer.sequence)
        return self._identity

    @property
    def contig_identity(self):
        """ Calculate and return the contiguous identity between the
        pair sequences. """
        if self._contig_identity is not None:
            return self._contig_identity

        matcher = SequenceMatcher()
        matcher.set_seq1(self.forward_primer.sequence)
        matcher.set_seq2(self.reverse_primer.sequence)
        _, _, size = matcher.find_longest_match(
            0, len(self.forward_primer.sequence),
            0, len(self.reverse_primer.sequence))
        self._contig_identity = int(size)
        return self._contig_identity

class Sequence:
    """ String wrapper to represent a DNA sequence. This is useful
    because it allows us to quickly generate complements, reverse
    complements, gc content percentages, and other characteristics. """

    def __init__(self, sequence):
        self.sequence = str(sequence)
        self._complementary_score = None
        self._contig_complementary_score = None
        self._gc = None
        self._tm = None
        

    def __eq__(self, other):
        if not isinstance(other, Sequence):
            return False
        return self.sequence == other.sequence

    def __getitem__(self, key):
        return Sequence(self.sequence.__getitem__(key))

    def __hash__(self):
        return hash(self.sequence)

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f'Sequence(sequence={self.sequence})'

    def complement(self):
        """ Return a Sequence object that is the complement of
        self.sequence. """
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return Sequence(''.join([comp[i] if i in comp else i for i in list(self.sequence)]))

    def reverse(self):
        """ Return a Sequence object that is the reverse of
        self.sequence. """
        return Sequence(self.sequence[::-1])

    def rev_comp(self):
        """ Return the reverse complement of self.sequence. """
        return self.complement().reverse()

    @property
    def gc(self):
        """ Returns a double in [0, 1] representing the ratio of G and
        C nucleotides. """
        if self._gc is not None:
            return self._gc

        if not self.sequence:
            return 0

        self._gc = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence)
        return self._gc

    @property
    def tm(self):
        """ Returns an estimated melting temperature of this sequence.
        Note: this is only accurate for sequences of ~16-24 characters.

        TODO: Include the sources for these values.
        """

        if self._tm is not None:
            return self._tm

        enthalpy = {'AA': -7.9, 'AT': -7.2, 'AC': -8.4, 'AG': -7.8,
                    'TA': -7.2, 'TT': -7.9, 'TC': -8.2, 'TG': -8.5,
                    'CA': -8.5, 'CT': -7.8, 'CC': -8.0, 'CG': -10.6,
                    'GA': -8.2, 'GT': -8.4, 'GC': -9.8, 'GG': -8.0}
        entropy = {'AA': -22.2, 'AT': -20.4, 'AC': -22.4, 'AG': -21,
                   'TA': -21.3, 'TT': -22.2, 'TC': -22.2, 'TG': -22.7,
                   'CA': -22.7, 'CT': -21.0, 'CC': -19.9, 'CG': -27.2,
                   'GA': -22.2, 'GT': -22.4, 'GC': -24.4, 'GG': -19.9}

        # Split the sequence into strings of length 2 and get their value.
        chunks = [self.sequence[i:i+2] for i in range(len(self.sequence)-1)]
        try:
            h = sum(enthalpy[chunk] for chunk in chunks)
            s = sum(entropy[chunk] for chunk in chunks)
            melting_temp = 1000*(h+2.4)/(s-0.7*len(self.sequence)+2-32.1)-273.15
            melting_temp = round(melting_temp, 2)
        except KeyError:
            melting_temp = float('nan')

        self._tm = melting_temp
        return self._tm

    @property
    def complementary_score(self):
        """ Return the self complementary score of this primer. """
        if isinstance(self._complementary_score, int):
            return self._complementary_score
            
        self._complementary_score = utils.complementary_score(self,
                                                              self.reverse())

        return self._complementary_score

    @property
    def contig_complementary_score(self):
        """ Complement the reverse sequence and check for equality. This is
        identical to checking for complementary of the reverse sequence.
        difflib.SequenceMatcher is very fast at this.

        Since this is an expensive function, the score is saved after the
        first calculation. Then, it is returned in any subsequent calls.
        """
        if self._contig_complementary_score is not None:
            return self._contig_complementary_score

        self._contig_complementary_score = utils.contig_complementary_score(self,
                                                                            self.reverse())

        return self._contig_complementary_score

    def has_contig_gc_at(self, gc, at):
        """ Return True if primer has gc contiguous G/C or at contiguous A/T. """
        pattern = re.compile('[GC]{' + str(gc) + '}|[AT]{' + str(at) + '}')
        if re.search(pattern, self.sequence):
            return True
        return False

    def has_in_last(self, gc, at, p):
        """
        Args:
            gc: The number of G/C nucleotides to search for.
            at: The number of A/T nucleotides to search for.
            p: The last p nucleotides are searched.

        Returns:
            True if this sequence has gc G/C or at A/T in the last p
            bases.

        Raises:
            IndexError: if p is greater than the sequence length.
        """
        if p > len(self.sequence):
            raise IndexError("Argument p is greater than the sequence length.")

        seq = self.sequence[-1*p:]

        A = seq.count('A')
        T = seq.count('T')
        C = seq.count('C')
        G = seq.count('G')

        return G+C >= gc or A+T >= at

    def has_mononucleotide_repeat(self, n):
        """
        Args:
            n: Desired length of mononucleotide repeat.

        Returns
            True if this sequence has a mononucleotide repeat.

        Raises:
            None
        """
        pattern = re.compile('A+|G+|C+|T+')
        result = re.findall(pattern, self.sequence)
        return len(max(result, key=len)) >= n

    def has_dinucleotide_repeat(self, n):
        """
        Args:
            n: Desired length of dinucleotide repeat.

        Returns:
            True if this sequence has a dinucleotide repeat of length n.

        Raises:
            None
        """
        n = str(n)
        pattern = re.compile('(AT){'+ n +'}|(TA){'+ n +'}|(AC){'+ n +'}|(CA){'+ n
                             +'}|(AG){'+ n +'}|(GA){'+ n +'}|(TG){'+ n +'}|(GT){'+ n
                             +'}|(TC){'+ n +'}|(CT){'+ n +'}|(GC){'+ n +'}|(CG){'
                             + n + '}')
        if re.search(pattern, self.sequence):
            return True
        return False

    def is_self_complementary(self, n, m):
        """
        Return true if the primer has contiguous complementary >= n or
        the difference between primer length and m <= m.

        Args:
            n: Length of desired contiguous complementary
            m: Number of non-complementary nucleotides.

        Returns:
            True: if the primer has contiguous complementary >= n or
                the difference between primer length and complementary
                score <= m.

        Raises:
            None
        """
        return (self.contig_complementary_score >= n
                or len(self) - self.complementary_score <= m)

class Primer(Sequence):
    """ Represents a PCR primer. 

    Args:
        sequence: The string of nucleotides that make up this primer.
        span: A tuple with start and end indices of the primer's binding
            site, where start is inclusive and end is exclusive.
        strand: 1 represents the plus strand, -1 is the minus strand.
        custom: Set to true if this is a custom primer the user
            specified.
    """

    def __init__(self, sequence, span, strand, custom=False):
        # All sequences are oriented 5'->3' on its corresponding strand.
        # But, start and end points are indexed based on the plus strand
        super().__init__(sequence)
        self.span = span
        self.strand = strand
        self.custom = custom

    def __repr__(self):
        return (f'Primer(sequence={str(self.sequence)}, '
                f'span={self.span}, '
                f'strand={self.strand}, '
                f'custom={self.custom})')

    @property
    def start(self):
        return self.span[0]

    @property
    def end(self):
        return self.span[1]
