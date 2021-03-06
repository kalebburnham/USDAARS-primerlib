import re
from copy import deepcopy

from .utils import (add_tails, complementary_score, contig_complementary_score,
                    cut, segregate, rtailed, rsorted)

class Sequence:
    """
    A sequence is a wrapper around a string of nucleotides that allows
    operations such as complements, reverse complements, and melting
    temperature estimations.

    Attributes:
        sequence: The string of nucleotides.
        gc: The GC content of this sequence. Between 0 and 1.
        tm: The melting temperature of this sequence, in Celsius.
        complementary_score: The number of complementary nucleotides
            this sequence has with its reverse.
        contig_complementary_score: The number of contigous
            complementary nucleotides this sequence has with its
            reverse.

    Public Methods:
        complement: Complement of this sequence.
        reverse: Reverse this sequence.
        rev_comp: Reverse complement this sequence.
        has_contig_gc_at: Check if this sequence has a contiguous G/C
            or A/T sequence.
        has_in_last: Check if this sequence has a specific number of
            G/C or A/T bases in the last X nucleotides.
        has_mononucleotide_repeat: Check if this sequence has a repeat
            of a single nucleotide.
        has_dinucleotide_repeat: Check if this sequence has a repeat of
            two nucleotides.
    """

    def __init__(self, sequence):
        self.sequence = str(sequence)

    def __add__(self, other):
        return Sequence('{}{}'.format(self.sequence, other))

    def __eq__(self, other):
        if not isinstance(other, Sequence):
            return False

        return self.sequence == other.sequence

    def __getitem__(self, key):
        return Sequence(self.sequence.__getitem__(key))

    def __hash__(self):
        return hash(self.sequence)

    def __iter__(self):
        for c in self.sequence:
            yield Sequence(c)

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
        return Sequence(''.join([comp[i] if i in comp else i for i in list(str(self.sequence))]))

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
        if not self.sequence:
            return 0
        return (str(self.sequence).count('G') + str(self.sequence).count('C')) / len(self.sequence)

    @property
    def tm(self):
        """Returns an estimated melting temperature of this sequence.
        Note: this is only accurate for sequences of ~16-24 characters.

        Source: docs/starp/Nearest-Neighbor Thermodynamic Parameters.pdf
        """
        if re.search('[^ACGT]', self.sequence) is not None:
            return float('nan') # Cannot calculate tm with bad characters.

        enthalpy = {'AA': -7.9, 'AT': -7.2, 'AC': -8.4, 'AG': -7.8,
                    'TA': -7.2, 'TT': -7.9, 'TC': -8.2, 'TG': -8.5,
                    'CA': -8.5, 'CT': -7.8, 'CC': -8.0, 'CG': -10.6,
                    'GA': -8.2, 'GT': -8.4, 'GC': -9.8, 'GG': -8.0}
        entropy = {'AA': -22.2, 'AT': -20.4, 'AC': -22.4, 'AG': -21,
                   'TA': -21.3, 'TT': -22.2, 'TC': -22.2, 'TG': -22.7,
                   'CA': -22.7, 'CT': -21.0, 'CC': -19.9, 'CG': -27.2,
                   'GA': -22.2, 'GT': -22.4, 'GC': -24.4, 'GG': -19.9}

        # Split the sequence into strings of length 2 and get their value.
        h = sum(list(enthalpy[self.sequence[i:i+2]]
                     for i in range(len(self.sequence)-1)))
        s = sum(list(entropy[self.sequence[i:i+2]]
                     for i in range(len(self.sequence)-1)))

        N = len(self.sequence) + 2
        tm = 1000*(h+2.4)/(s-0.7*N-32.1)-273.15
        return tm

    @property
    def complementary_score(self):
        """ Return the self complementary score of this sequence. """
        if hasattr(self, '_complementary_score'):
            return self._complementary_score

        setattr(self, '_complementary_score',
                complementary_score(self, self.reverse()))

        return self._complementary_score

    @property
    def contig_complementary_score(self):
        """ Return the self contiguous complementary score of this
        primer. """
        if hasattr(self, '_contig_complementary_score'):
            return self._contig_complementary_score

        setattr(self, '_contig_complementary_score',
                contig_complementary_score(self, self.reverse()))

        return self._contig_complementary_score

    def has_repeated_nucleotide(self, n):
        """ Returns True if primer has >= n of a single nucleotide. """

        seq = self.sequence
        num_A = seq.count('A')
        num_T = seq.count('T')
        num_C = seq.count('C')
        num_G = seq.count('G')

        return num_A >= n or num_T >= n or num_C >= n or num_G >= n

    def has_contig_gc_at(self, num_gc, num_at):
        """ Return True if primer has num_gc contiguous G/C or num_at
        contiguous A/T. """
        pattern = re.compile('[GC]{' + str(num_gc) + '}|[AT]{' + str(num_at) + '}')
        return bool(re.search(pattern, self.sequence))

    def has_in_first(self, num_gc, num_at, p):
        """
        Args:
            num_gc: The number of G/C nucleotides to search for.
            num_at: The number of A/T nucleotides to search for.
            p: The first p nucleotides are searched.

        Returns:
            True if this sequence has gc G/C or at A/T in the first p
            bases.

        Raises:
            IndexError: if p is greater than the sequence length.
        """
        if p > len(self.sequence):
            raise IndexError("Argument p is greater than the sequence length.")

        seq = self.sequence[:p]

        A = seq.count('A')
        T = seq.count('T')
        C = seq.count('C')
        G = seq.count('G')

        return G+C >= num_gc or A+T >= num_at

    def has_in_last(self, num_gc, num_at, p):
        """
        Args:
            num_gc: The number of G/C nucleotides to search for.
            num_at: The number of A/T nucleotides to search for.
            p: The last p nucleotides are searched.

        Returns:
            True if this sequence has gc G/C or at A/T in the last p
            bases. Else, False.

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

        return G+C >= num_gc or A+T >= num_at

    def has_mononucleotide_repeat(self, num_gc, num_at):
        """
        Args:
            num_gc: Desired length of G or C repeat.
            num_at: Desired length of A or T repeat.

        Returns
            True if this sequence has a mononucleotide repeat. Else, False.
        """
        num_gc = str(num_gc)
        num_at = str(num_at)
        pattern = re.compile('G{' + num_gc + '}|C{' + num_gc + '}|A{'
                             + num_at + '}|T{' + num_at + '}')
        return bool(re.search(pattern, self.sequence))

    def has_dinucleotide_repeat(self, n):
        """
        Args:
            n: Desired length of dinucleotide repeat.

        Returns:
            True if this sequence has a dinucleotide repeat of length n.
        """
        n = str(n)
        pattern = re.compile('(AT){'+ n +'}|(TA){'+ n +'}|(AC){'+ n +'}|(CA){'+ n
                             +'}|(AG){'+ n +'}|(GA){'+ n +'}|(TG){'+ n +'}|(GT){'+ n
                             +'}|(TC){'+ n +'}|(CT){'+ n +'}|(GC){'+ n +'}|(CG){'
                             + n + '}')
        return bool(re.search(pattern, self.sequence))

    @staticmethod
    def hamming(s1, s2):
        """Return the Hamming distance between equal-length sequences.
        Source: Wikipedia """
        s1 = str(s1)
        s2 = str(s2)
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length.")
        return sum(el1 != el2 for el1, el2 in zip(s1, s2))

class Primer(Sequence):
    """ Represents a PCR primer.
    
    Attributes:
        sequence: Inherited from :class:Sequence.
        allele1_span: The tuple [start, end) of this primer's binding site on
            the first allele.
        allele2_span: The tuple [start, end) of this primer's binding site on
            the second allele.
        strand: Either 1 or -1 specifying if this primer resides on the plus
            or minus strand.
        allele1_start: The starting position of this Primer on allele1.
        allele1_end: The 'end' position of this Primer on allele1 from the span.
        allele2_start: The starting position of this Primer on allele2.
        allele2_end: The 'end' position of this Primer on allele2 from the span.

    Public Methods:
        rev_comp: Reverse complement this primer while switching strands.
    """

    def __init__(self, sequence: str, allele1_span: tuple, allele2_span: tuple, strand):
        # All sequences are oriented 5'->3' on its corresponding strand.
        # But, start and end points are indexed based on the plus strand
        super().__init__(sequence)
        self.allele1_span = allele1_span
        self.allele2_span = allele2_span
        # Strand 1 corresponds to the plus strand. Strand -1 is the minus strand.
        self.strand = strand

    def __repr__(self):
        return (f'Primer(sequence={str(self.sequence)}, '
                f'allele1_span={self.allele1_span}, '
                f'allele2_span={self.allele2_span}, '
                f'strand={self.strand})')

    def __len__(self):
        return len(self.sequence)

    def __contains__(self, item):
        return item in self.sequence

    def __eq__(self, other):
        return (self.sequence == other.sequence
                and self.strand == other.strand
                and self.allele1_span == other.allele1_span
                and self.allele2_span == other.allele2_span)

    def __getitem__(self, key):
        return self.sequence[key]

    def __hash__(self):
        return hash((self.sequence, self.allele1_span, self.allele2_span, self.strand))

    @property
    def allele1_start(self):
        return self.allele1_span[0]

    @property
    def allele1_end(self):
        return self.allele1_span[1]

    @property
    def allele2_start(self):
        return self.allele2_span[0]

    @property
    def allele2_end(self):
        return self.allele2_span[1]

    def rev_comp(self):
        """ Reverse complement sequence and return a new primer on the
        opposite strand. """
        seq = Sequence(self.sequence)
        return Primer(seq.rev_comp(), allele1_span=self.allele1_span,
                      allele2_span=self.allele2_span, strand=self.strand*-1)
    
class AmasPrimer(Sequence):
    """
    A container for an AMAS primer.

    Attributes:
        tail: The added tail on the 5' end.
        allele_num: The allele this primer originated from. The first
            allele entered is allele 1 and the second is allele 2.
        span: The beginning and ending positions of this primer. Should
        be [inclusive, exclusive)
        start: The first number in the span.
        end: The second number in the span.

    Public Methods:
        tm: The melting temperature of this primer's original sequence
            (before substitution) without the tail.
        rev_comp: Reverse complement everything except the tail. Tails
            are never rev_comp'ed.
        html: Return a string of markup with the tail underlined and
            italicized. This is used in the HTML templates.
    """

    def __init__(self, sequence, allele_num, span, strand=1, original_seq=None):
        super().__init__(sequence)
        self.tail = Sequence('')  # A sequence object.
        self.allele_num = None  # 1 or 2
        self.span = span
        self.strand = strand

        # The Tm value is based on the original strand. Substitution
        # modifies the primer sequence so the original is kept track
        # of here.
        if original_seq:
            self._original_seq = original_seq
        else:
            self._original_seq = self.sequence

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return str(self.tail) + str(self.sequence)

    def __repr__(self):
        return (f'AmasPrimer('
                f'tail={str(self.tail)}, '
                f'sequence={str(self.sequence)}, '
                f'allele_num={self.allele_num}, '
                f'span={self.span}, '
                f'strand={self.strand}'
                f')')

    @property
    def start(self):
        return self.span[0]
    
    @property
    def end(self):
        return self.span[1]

    @property
    def tm(self):
        """ Return the melting temperature for the original sequence
        without the tail. """
        return Sequence(self._original_seq).tm

    def rev_comp(self):
        """ Return a new AMAS primer with only the sequence changed.
        Tails are not reverse complemented."""
        seq = Sequence(self.sequence).rev_comp()
        original_seq = Sequence(self._original_seq).rev_comp()
        primer = AmasPrimer(seq, self.allele_num, self.span, strand=self.strand*-1, original_seq=original_seq)
        primer.tail = self.tail
        return primer
    
    def html(self):
        """ One document shows the tail underlined and the substitutions
        highlighted. """
        markup = f'<span><u><i>{str(self.tail)}</i></u>{str(self.sequence)}</span>'
        return markup

class Snp:
    """
    An object representing a single-nucleotide polymorphism.

    The only acceptable SNPs are substitutions, insertions, and
    deletions.

    Substitution: '{prefix}.{index}{ref_nucleotide}>{new_nucleotide}
        Ex. 'SequenceA.16C>G'
    Insertion: '{prefix}.{index}ins{new_nucleotide}
        Ex. 'SequenceB.178insA
    Deletion: '{prefix}.{index}del
        Ex. 'SequenceC.19del'

    Attributes:
        descriptor: A standardized string describing the SNP. See
            http://varnomen.hgvs.org/recommendations/DNA/.
        position: The zero-indexed position of the SNP.
        type: A string saying the type of the SNP. Current possible
            values are "substitution", "insertion", and "deletion".
        ref_nucleotide: The first possible nucleotide in the SNP. For
            example, in the SNP [C/G], ref_nucleotide is C.
        new_nucleotide: The second possible nucleotide in the SNP. For
            example, in the SNP [C/G], new_nucleotide is G.
        nucleotides: A set containing the two nucleotides.
            Usage:
            >> snp.ref_nucleotide == 'G'
            >> snp.new_nucleotide == 'C'
            >> snp.nucleotides == {C, G}
            True

    """

    def __init__(self, descriptor):
        """ It is assumed the descriptor has zero-indexed positions.
        Currently this object should accept the Human Genome Variation
        Society SNP standard as described at
        http://varnomen.hgvs.org/recommendations/DNA/.
        The only accepted SNPs currently are substitutions, deletions,
        and insertions.

        Descriptor Examples:
            Substitution: '.100A>C'
            Insertion: 'SequenceA.87insG'
            Deletion: 'SequenceB.152del'

        Note that the descriptor for deletions contains no information
        on what the deleted nucleotide is. This information must be
        entered manually.
        """
        self.descriptor = descriptor

        parser = SnpParser(self.descriptor)
        self.prefix = parser.prefix()
        self.position = parser.position()
        self.type = parser.type()
        self.ref_nucleotide = parser.ref_nucleotide()
        self.new_nucleotide = parser.new_nucleotide()

    def __eq__(self, other):
        return (self.descriptor == other.descriptor
                and self.prefix == other.prefix
                and self.position == other.position
                and self.type == other.type
                and self.ref_nucleotide == other.ref_nucleotide
                and self.new_nucleotide == other.new_nucleotide)

    def __hash__(self):
        return super.__hash__(self)

    def __str__(self):
        return self.descriptor

    def __repr__(self):
        return (f'Snp(descriptor={self.descriptor}, '
                f'type={self.type}, '
                f'prefix={self.prefix}, '
                f'position={self.position}, '
                f'ref_nucleotide={self.ref_nucleotide}, '
                f'new_nucleotide={self.new_nucleotide})')

    @property
    def nucleotides(self):
        # Type frozen set so it is hashable
        return frozenset((self.ref_nucleotide, self.new_nucleotide))

class SnpParser:
    """
    Given a descriptor, this class parses out the relevant information
    for quicker access.
    """
    def __init__(self, descriptor):
        self.descriptor = descriptor

    def position(self):
        """
        Extract the start index from self.descriptor.

        Args:
            None

        Returns:
            An int representing the start index of self.descriptor.
        """
        suffix = self.descriptor.replace(self.prefix(), '')
        match = re.match(r'\d+', suffix)
        return int(match.group())

    def type(self):
        """ Current supported types are substitutions, insertions,
        and deletions.
        """
        suffix = self.descriptor.replace(self.prefix(), '')
        if re.search('>', suffix):
            snp_type = 'substitution'
        elif re.search('ins', suffix):
            snp_type = 'insertion'
        elif re.search('del', suffix):
            snp_type = 'deletion'
        else:
            snp_type = ''

        return snp_type

    def prefix(self):
        """ Match everything up to the colon (if it exists),
        then continue matching up to and including the next period.
        Example: NC_000023.10:g.32867861_32867862insT
        Prefix is 'NC_000023.10:g.'.
        """
        pattern = r'[^:]*[:]?[A-Za-z\d]*\.'
        return re.match(pattern, self.descriptor).group()

    def ref_nucleotide(self):
        """ If the type isn't a substitution, we don't know what the
        original nucleotide is.
        Example: g.123A>G - ref nucleotide is A """

        suffix = self.descriptor.replace(self.prefix(), '')
        if self.type() == 'substitution':
            # Match the letter right before the '>'
            pattern = re.compile('[GCAT]{1}(?=>)')
            nucleotide = re.search(pattern, suffix).group()
        else:
            nucleotide = ''

        return nucleotide

    def new_nucleotide(self):
        """ Example: g.123A>G - new nucleotide is G """
        suffix = self.descriptor.replace(self.prefix(), '')
        if self.type() == 'substitution':
            # Match the letter following the '>'
            pattern = re.compile('(?<=>)[GCAT]{1}')
            nucleotide = re.search(pattern, suffix).group()
        elif self.type() == 'insertion':
            pattern = re.compile('(?<=ins)[GCAT]+')
            nucleotide = re.search(pattern, suffix).group()
        else:
            nucleotide = ''

        return nucleotide

class StarpGroup:

    """
    A class that contains primer and snp information for a single
    location. 

    Attributes:
        amas1: The AMAS primer on the first allele.
        amas2: The AMAS primer on the second allele.
        snp_position: The focus snp is either at the 'first' or 'last'
            nucleotide in the amas primers. If the snp position is
            'first', rprimers will be looked for before the amas
            primers. If it is last, they will be searched for after the
            amas primers.
        rcandidates: All possible reverse primers. From this list, the
            best primers are chosen to pair with the AMAS primers.
            These primers should have already been sorted, filtered by
            their contents, and filtered by their binding sites.
        rprimers: The primers that are considered 'paired' with the AMAS
            primers. These are the primers that are displayed.
        snp: The Snp object which the primers are being constructed about.
        num_rprimers: The maximum number of primers to pair with the
            AMAS primers. Defaults to 3.
    """

    def __init__(self, amas1, amas2, snp_position, rcandidates=[], snp=None, num_rprimers=3):
        self.amas1 = amas1
        self.amas2 = amas2
        self.snp_position = snp_position
        self.rcandidates = rcandidates
        self.rprimers = []
        self.snp = snp
        self.num_rprimers = num_rprimers

    def substitute_amas_bases(self):
        """ Substitutes AMAS primer bases according to Long and
        reassigns amas1 and amas2. 

        Returns:
            [Modified amas1, Modified amas2]
        """
        from .amasfactory import (substitute_bases,)

        self.amas1.sequence, self.amas2.sequence = substitute_bases([self.amas1.sequence, self.amas2.sequence], snp_position=self.snp_position)
        return [self.amas1, self.amas2]
    
    def set_rprimers(self):
        """ Sorts the rcandidates and selects the first 'num_rprimers' primers
        that are in a valid position according to the AMAS primers and will
        not bind to them.

        Returns:
            The best rprimers matching this pair of AMAS primers.
        """
        sorted_primers = rsorted(self.rcandidates)

        for primer in sorted_primers:

            # If this primer is not in a valid position or it overlaps with an
            # AMAS primer, move on to the next one.
            if self.snp_position == 'first':
                if primer.allele1_end > self.amas1.start:
                    continue
                if primer.allele2_end > self.amas2.start:
                    continue
            else:
                if primer.allele1_start < self.amas1.end:
                    continue
                if primer.allele2_start < self.amas2.end:
                    continue

            # Check complementarity with AMAS primers. If this r primer
            # has at least 5 unpaired nucleotides with both AMAS primers,
            # add it to the list.
            if (len(primer) - complementary_score(primer, self.amas1) > 5
                    and len(primer) - complementary_score(primer, self.amas2) > 5):
                self.rprimers.append(primer)

            if len(self.rprimers) >= self.num_rprimers:
                break

        return self.rprimers

    def segregate(self, tail1, tail2):
        """ Splits this group into max two other groups. In one, amas1
        and amas2 have tail1 and tail2 respectively. The corresponding
        rprimers agree with these tails. By 'agree', it is meant that
        the tails added to the AMAS primers are the same tails that are
        defined by the rprimers according to
        'docs/starp/STARP R primer design[4311].docx', page 6.

        In the other group, amas1 and amas2 have tail2 and tail1
        respectively. Again, the group's rprimers agree with these
        tails.

        Once the groups has been made, the tails need to be cut if they
        have overlapping bases at the 3' end.
        """
        
        group1 = StarpGroup(deepcopy(self.amas1), deepcopy(self.amas2), self.snp_position, snp=self.snp)
        group1.amas1.tail = tail1
        group1.amas2.tail = tail2

        # Choose the rprimers that lead to tail1 attached to amas1 and tail2
        # attached to amas2.
        group1.rprimers = list(filter(lambda rprimer: self.assigns_tail(
                                             group1.amas1.tail, group1.amas2.tail,
                                             group1.amas1, group1.amas2, rprimer, self.snp_position),
                                      self.rprimers))


        group2 = StarpGroup(deepcopy(self.amas1), deepcopy(self.amas2), self.snp_position, snp=self.snp)
        group2.amas1.tail = tail2
        group2.amas2.tail = tail1
        # Choose the rprimers that lead to tail2 attached to amas1 and tail1
        # attached to amas2.
        group2.rprimers = list(filter(lambda rprimer: self.assigns_tail(
                                             group2.amas1.tail, group2.amas2.tail,
                                             group2.amas1, group2.amas2, rprimer, self.snp_position),
                                      self.rprimers))

        # Cut the ends of the tails if they have overlapping bases with the
        # AMAS primers as instructed by
        # 'docs/starp/STARP R primer design[4311].docx', page 6.
        group1.amas1.tail = cut(group1.amas1.tail, Sequence(group1.amas1.sequence), snp_position=group1.snp_position)
        group1.amas2.tail = cut(group1.amas2.tail, Sequence(group1.amas2.sequence), snp_position=group1.snp_position)

        group2.amas1.tail = cut(group2.amas1.tail, Sequence(group2.amas1.sequence), snp_position=group2.snp_position)
        group2.amas2.tail = cut(group2.amas2.tail, Sequence(group2.amas2.sequence), snp_position=group2.snp_position)

        # Note: none of these groups may not have any rprimers.
        return [group1, group2]

    def assigns_tail(self, amas1_tail, amas2_tail, amas1, amas2, rprimer, snp_position):
        """ Returns true if the tails defined by the rprimer are equal to
        amas1_tail and amas2_tail.

        Source: docs/starp/STARP R primer design[4311].docx page 6.

        Args:
            amas1_tail: The tail attached to amas1.
            amas2_tail: The tail attached to amas2.
            amas1: The AmasPrimer object on the first allele.
            amas2: The AmasPrimer object on the second allele.
            rprimer: The Primer object that is being checked.
            snp_position: A string specifying the relative position of
                the SNP between these AMAS primers. Either 'first' or
                'last'.
        """
        tail1 = Sequence('GCAACAGGAACCAGCTATGAC')
        tail2 = Sequence('GACGCAAGTGAGCAGTATGAC')

        amplicon1 = abs(amas1.start - rprimer.allele1_end)
        amplicon2 = abs(amas2.start - rprimer.allele2_end)

        if snp_position == 'first':
            if amas1.strand == 1:
                amas1 = amas1.rev_comp()
            if amas2.strand == 1:
                amas2 = amas2.rev_comp()

        return_val = False

        # The tail defined by the rprimer depends on the amplicon lengths
        # between it and the AMAS primers. See
        # docs/starp/STARP R primer design[4311].docx page 6
        if amplicon1 - amplicon2 >= 8:
            return_val = (tail1 == amas1_tail and amas2 == amas2_tail)
        elif amplicon1 - amplicon2 >= 1:
            return_val = (tail2 == amas1_tail and tail1, amas2 == amas2_tail)
        elif amplicon1 - amplicon2 == 0:
            # This branch comes from Table 3 in the STARP paper.
            assigned_tail = {
                frozenset({'C', 'G'}): {'C': 1, 'G': 2},
                frozenset({'C', 'T'}): {'C': 1, 'T': 2},
                frozenset({'C', 'A'}): {'C': 1, 'A': 2},
                frozenset({'G', 'T'}): {'G': 1, 'T': 2},
                frozenset({'G', 'A'}): {'G': 1, 'A': 2},
                frozenset({'T', 'A'}): {'T': 1, 'A': 2}
            }

            nucleotides = frozenset({str(amas1[-1]), str(amas2[-1])})

            if assigned_tail[nucleotides][str(amas1[-1])] == 1:
                return_val = (tail1 == amas1_tail and tail2 == amas2_tail)
            else:
                return_val = (tail2 == amas1_tail and tail1 == amas2_tail)

        elif amplicon1 - amplicon2 >= -7:
            return_val = (tail1 == amas1_tail and tail2 == amas2_tail)
        else:
            return_val = (tail2 == amas1_tail and tail1, amas2 == amas2_tail)

        return return_val
