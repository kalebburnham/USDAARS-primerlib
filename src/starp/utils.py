"""
License information goes here.
"""

from difflib import SequenceMatcher
import itertools

from Bio import Align
import regex

from .exceptions import StarpError

def add_tails(amas1, amas2, amplicon1, amplicon2, snp_position):
    """
    Add tails to AMAS primers according to 'STARP R primer design[4311].docx'.
    Note that downstream primers must be reverse complemented before
    adding the unmodified tail to the 5' end.

    Args:
        amas1: An AmasPrimer object.
        amas2: An AmasPrimer object.
        amplicon1: The length of the amplicon from the first allele.
        amplicon2: The length of the amplicon from the second allele.
        snp_position: Position of the SNP. Should be either 'first'
            or 'last'.
    """
    from .models import Sequence

    tail1 = Sequence('GCAACAGGAACCAGCTATGAC')
    tail2 = Sequence('GACGCAAGTGAGCAGTATGAC')

    if snp_position == 'first':
        if amas1.strand == 1:
            amas1 = amas1.rev_comp()
        if amas2.strand == 1:
            amas2 = amas2.rev_comp()

    if amplicon1 - amplicon2 >= 8:
        amas1.tail = cut(tail1, amas1)
        amas2.tail = cut(tail2, amas2)
    elif amplicon1 - amplicon2 >= 1:
        amas1.tail = cut(tail2, amas1)
        amas2.tail = cut(tail1, amas2)
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
            amas1.tail = cut(tail1, amas1)
            amas2.tail = cut(tail2, amas2)
        else:
            amas1.tail = cut(tail2, amas1)
            amas2.tail = cut(tail1, amas2)

    elif amplicon1 - amplicon2 >= -7:
        amas1.tail = cut(tail1, amas1)
        amas2.tail = cut(tail2, amas2)
    else:
        amas1.tail = cut(tail2, amas1)
        amas2.tail = cut(tail1, amas2)

    return (amas1, amas2)

def add_rtails(rprimers):
        """ Add tails to rprimers with low melting temperature. """
        low, high = segregate(rprimers)
        low = rtailed(low)
        return low + high

def aligned_position(pos, snps):
    """
    Returns an updated position after the dashes are added back into
    the sequence.
    Snp positions are given relative to allele1 without any dashes,
    Example:
        With snp = '.2insA', allele1 = 'ATGCACTG', and pos=4, the
        aligned sequence is 'AT-GCACTG' and the new position is 5
        because we added a single dash before it.
    """

    insertion_count = len([snp for snp in snps
                           if snp.position < pos and snp.type == 'deletion'])

    return pos + insertion_count

def cut(tail, primer):
    """
    Checks the overlap between the 3' end of 'TATGAC' and the 5' end of
    the primer, and returns the tail with the overlapping bases
    removed.

    For example, with tail='GCAACAGGAACCAGCTATGAC' and
    primer='TGACGGATCTAGACTACTGACGTCA', 'TGAC' overlaps them. So,
    the :class:Sequence 'GCAACAGGAACCAGCTA' is returned.

    Args:
        tail: Should be one of the AMAS tails of :class:Sequence.
            These are 'GCAACAGGAACCAGCTATGAC' and
            'GACGCAAGTGAGCAGTATGAC'.
        primer: The :class:Primer primer to check overlaps on its 5'
            end.
    """
    s = 'TATGAC'

    overlap = 0
    for idx in range(6):
        if s[5-idx:] == primer.sequence[:idx]:
            overlap = idx+1

    if overlap == 0:
        return tail
    else:
        return tail[:0-overlap]

def binding_sites(sequences: tuple, primer, stop=2):
    """
    Returns a list of potential binding sites of this primer on these
    sequences. A binding site is defined as a match with <= 4 nucleotide
    subsitutions with a mismatch on the 3' end AND 2 more mismatches on
    the 2nd, 3rd, or 4th position from the 3' end.

    Note that this method does not specify which sequence a binding site
    was found on, only its position within a sequence. It is usually
    only important to know if a primer has more than one binding site,
    not necessarily where it is located. In the case of custom primers,
    it will be assumed that the user will give a primer on the reference
    sequence. Thus, we need to know the location of its binding site.
    This is given by the regex.Match object.

    Since this is an expensive operation, the 'stop' argument is included
    to immediately exit if the specified number of binding sites are
    found. This saves time and resources when dealing with repetitive
    sequences.

    Args:
        sequences: Reference sequences and high scoring pair sequences
            to check for binding sites.
        primer: Look for binding sites of this primer.
        stop: Stop looking for binding sites if this many have already
            been found.

    Returns:
        A list of regex.Match objects representing binding sites for the
        primer.

    Raises:
        None
    """
    TOLERANCE = 4  # Number of acceptable nucleotide substitutions.

    sites = list()
    pattern = '(' + primer.sequence + '){s<=' + str(TOLERANCE) + '}'
    revcomp_pattern = '(' + primer.rev_comp().sequence + '){s<=' + str(TOLERANCE) + '}'

    for sequence in sequences:
        matches_iter = regex.finditer(pattern, str(sequence),
                                      overlapped=True)

        # This code is duplicated here because regex.finditer returns a
        # _regex.Scanner which is causing problems with itertools.chain
        # Namely, a StopIteration leaks out of the iterator so Python
        # complains that an object (a StopIteration) of type 'type' is
        # not iterable.
        for match in matches_iter:
            matched_seq = match.group()
            if (primer.sequence[-1] == matched_seq[-1]
                    or hamming(primer.sequence[-4:-1], matched_seq[-4:-1]) < 2):
                sites.append(match)

            if len(sites) >= stop:
                return sites

        matches_iter = regex.finditer(revcomp_pattern, str(sequence),
                                      overlapped=True)
        for match in matches_iter:
            if len(sites) >= stop:
                return sites

            matched_seq = match.group()
            if (primer.sequence[-1] == matched_seq[-1]
                    or hamming(primer.rev_comp().sequence[-4:-1], matched_seq[-4:-1]) < 2):
                sites.append(match)

        matches_iter = itertools.chain

    return sites

def complementary_score(s1, s2) -> int:
    """ Using Bio.Align from BioPython, calculate the score of a local
    sequence alignment using a very low open gap score as to completely
    eliminate gaps. Complement the second sequence because it is easier
    to check for same-ness than complementary nucleotides.

    ** Sequences should be oriented the same direction.

    Args:
        s1, s2: Strings of nucleotides.

    Returns:
        The maximum number of complementary nucleotides between s1 and
        s2 in a maximal alignment.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -1000 # Don't want any gaps!
    # score = the number of complementary nucleotides between s1 and s2.
    return aligner.score(str(s1), str(s2.complement()))

def contig_complementary_score(s1, s2) -> int:
    """ Complement the second sequence and check for equality. This is
    identical to checking for complementary of the reverse sequence.
    difflib.SequenceMatcher is very fast at this.

    ** Sequences should be oriented the same direction.

    Args:
        s1, s2: Sequences of nucleotides.

    Returns:
        The maximum number of contiguous complementary nucleotides
        between s1 and s2 in a maximal alignment.
    """
    matcher = SequenceMatcher()
    matcher.set_seq1(str(s1))
    matcher.set_seq2(str(s2.complement()))
    _, _, size = matcher.find_longest_match(0, len(s1), 0, len(s2))
    return size

def rgenerate(allele1, allele2, min_length, max_length):
    """
    Return a list of all possible R primers common to allele1 and
    allele2. However, the primers may not have a unique binding site on
    each allele.

    This algorithm creates a window of size 'size' and traverses
    allele1 and allele2 at the same time, while the span of each
    primer in allele1 is maintained in allele1_span and the span in
    allele2 is maintained in allele2_span. Stricly speaking, the spans
    are not accurate when there is a '-' in the sequence, but get
    updated when a dash is the first character in each sequence.

    Since a SNP cannot look like (-/-), it is not possible for two
    sequences to be equal if either contains a dash.

    Args:
        allele1: The first aligned allele.
        allele2: The second aligned allele.
        min_length: Minimum length of primer.
        max_length: Maximum length (inclusive) of primer.

    Returns:
        A list of all possible r primers common to allele1 and allele2.
    """
    from .models import Primer

    if len(allele1) != len(allele2):
        raise ValueError('Aligned alleles must be the same length.')

    allele1 = str(allele1)
    allele2 = str(allele2)
    candidates = []

    for size in range(min_length, max_length+1):
        allele1_span = (0, size)
        allele2_span = (0, size)

        for i in range(len(allele1)-size):
            if allele1[i:i+size] == allele2[i:i+size]:
                candidates.append(Primer(sequence=allele1[i:i+size],
                                         allele1_span=allele1_span,
                                         allele2_span=allele2_span,
                                         strand=1))

            if allele1[i] != '-':
                allele1_span = (allele1_span[0]+1, allele1_span[1]+1)

            if allele2[i] != '-':
                allele2_span = (allele2_span[0]+1, allele2_span[1]+1)

    return candidates

def rfilter_by_binding_sites(r_primers, allele1, allele2, nontargets):
    """
    allele1: The aligned first allele.
    allele2: The aligned second allele.

    TODO: Complete documentation.
    """

    candidates = []

    for primer in r_primers:
        # The primer should have at most 1 binding site in each allele.
        if (len(binding_sites((allele1,), primer)) > 1
                or len(binding_sites((allele2,), primer)) > 1):
            continue

        # Move on if the primer has other potential binding sites in
        # these sequences.

        if binding_sites(nontargets, primer, stop=1):
            continue

        if binding_sites(nontargets, primer.rev_comp(), stop=1):
            continue

        # The primer's binding sites are acceptable.
        candidates.append(primer)

    return candidates

def rtailed(rprimers: list) -> list:
    """
    Add tails to reverse primers as described in
    docs/Starp R primer design[4311].

    Then, returns the primers shorter than 28 bases and those with
    tm <= 62 degrees C.

    Args:
        r_primers: The primer list to add tails.

    Returns:
        A list of tailed R primers with length < 28 bp and tm <= 62
    """
    from .models import Primer

    tails = {'', 'C', 'G', 'CG', 'GC', 'CGC', 'GCG'}

    tailed = list()  # Primers with tails.
    for tail, primer in itertools.product(tails, rprimers):
        # Should the spans be updated when a tail is added?
        if primer.strand == 1:
            primer.sequence = primer.sequence + tail
        elif primer.strand == -1:
            if (primer.allele1_start - len(tail) >= 0
                    and primer.allele2_start - len(tail) >= 0):
                primer.sequence = tail + primer.sequence

    return [primer for primer in rprimers
            if len(primer) < 28
            and primer.tm <= 62]

def rfilter(r_primers: list) -> list:
    """
    Removes the R primers meeting any of the following conditions:
    - Contains invalid characters. The only accepted alphabet
      is {A, C, G, T}.
    - Has >= 10 contiguous G/C or >= 12 contiguous A/T
    - Has >= 8 As, Ts, Gs, or Cs
    - Has >= 6 di-nucleotide repeats
    - GC content is > 80% or < 20%

    Args:
        primers: Reverse primers. Strand doesn't matter.

    Returns:
        The primers that do not meet the above conditions.
    """

    candidates = [primer for primer in r_primers
                  if (not regex.search('[^ACGT]', primer.sequence)
                      and not primer.has_contig_gc_at(10, 12)
                      and primer.sequence.count('A') < 8
                      and primer.sequence.count('T') < 8
                      and primer.sequence.count('G') < 8
                      and primer.sequence.count('C') < 8
                      and not primer.has_dinucleotide_repeat(6)
                      and 0.20 <= primer.gc <= 0.80)]

    return candidates

def rfilter_tailed_primers(rprimers: list) -> list:
    """
    Removes the R primers meeting any of the following conditions:
    - Contains 28+ nucleotides bases
    - Has melting temperature > 62 degrees Celsius.
    - Contains 10+ contiguous self-complementary nucleotides.
    - Primer length - self complementary <= 4
    """

    candidates = [primer for primer in rprimers
                  if (len(primer) < 28
                      and primer.tm <= 62
                      and primer.contig_complementary_score < 10
                      and len(primer) - primer.complementary_score > 5)]

    return candidates

def segregate(primers: list):
    """
    Splits the primers into two groups. Group 1 contains primers
    with melting temperature between within [53 C, 58 C) and Group
    2 contains those with melting temperature within [58 C, 62 C].

    Args:
        primers: A list of primers to segregate.

    Returns:
        Two lists representing the two groups of primers stated
        above.
    """

    group1 = list(filter(lambda primer: 53 <= primer.tm < 58, primers))
    group2 = list(filter(lambda primer: 58 <= primer.tm <= 62, primers))

    return group1, group2

def rsorted(primers: list) -> list:
    """
    Sorts the primers with the following conditions:

    1) Has 9 contiguous (G and/or C) or 11 contiguous (A and/or T)
    2) Has 7 As, Ts, Gs, or Cs
    3) Has 5 di-nucleotide (AG, AC, TG, TC, GA, GT, CA, CT) Repeats
    4) Has (GC% > 75% or GC% < 25%)
    5) Has ≥ 8 contiguous complementarity or (primer length - max complementarity) ≤ 6
    6) Has 8 contiguous (G and/or C) or ≥ 9 contiguous (A and/or T)
    7) Has 6 As, Ts, Gs, or Cs
    8) Has 4 di-nucleotide (AG, AC, TG, TC, GA, GT, CA, CT) Repeats
    9) Has (GC% > 70% or GC% < 30%)
    10) Has 7 contiguous complementarity or (primer length - max complementarity) ≤ 8
    11) Has ≥ 6 contiguous (G and/or C) or ≥ 7 contiguous (A and/or T)
    12) Has 5 As, Ts, Gs, or Cs
    13) Has 3 di-nucleotide (AG, AC, TG, TC, GA, GT, CA, CT) Repeats
    14) Has (GC% > 65% or GC% < 35%)
    15) Has 6 contiguous complementarity or (primer length - max complementarity) ≤ 10
    16) Has 6 A/T or 5 G/C in the last seven bases
    17) Has 4 A/T or 3 G/C in the last four bases
    18) Has 5 contiguous complementarity or (primer length - max complementarity) ≤ 12
    19) Has ≥ 4 contiguous (G and/or C) or ≥ 5 contiguous (A and/or T)
    20) Has 4 As, Ts, Gs, or Cs
    21) Has (GC% > 60% or GC% < 40%)
    22) Has 4 contiguous complementarity or (primer length - max complementarity) ≤ 14
    23) Closest to 22 nucleotides

    Args:
        The R primers to be sorted.

    Returns:
        A sorted list of R primers with the best primers at the front.
    """

    # True values go towards the end of the list, so the best primers
    # are kept in the front.
    return sorted(primers,
                  key=lambda primer: (primer.has_contig_gc_at(9, 11),  # 1
                                      primer.has_repeated_nucleotide(7),  # 2
                                      primer.has_dinucleotide_repeat(5),  # 3
                                      primer.gc < 0.25 or primer.gc > 0.75,  # 4
                                      primer.contig_complementary_score >= 8,  # 5a
                                      len(primer) - primer.complementary_score <= 6,  #5b
                                      primer.has_contig_gc_at(8, 9),  # 6
                                      primer.has_repeated_nucleotide(6),  # 7
                                      primer.has_dinucleotide_repeat(4),  # 8
                                      primer.gc < 0.30 or primer.gc > 0.70,  # 9
                                      primer.contig_complementary_score >= 7,  # 10a
                                      len(primer) - primer.complementary_score <= 8,  # 10b
                                      primer.has_contig_gc_at(6, 7),  # 11
                                      primer.has_repeated_nucleotide(5),  # 12
                                      primer.has_dinucleotide_repeat(3),  # 13
                                      primer.gc < 0.35 or primer.gc > 0.65,  # 14
                                      primer.contig_complementary_score >= 6,  # 15a
                                      len(primer) - primer.complementary_score <= 10,  # 15b
                                      primer.has_in_first(5, 6, 7),  # 16
                                      primer.has_in_first(3, 4, 4),  # 17
                                      primer.contig_complementary_score >= 5,  # 18a
                                      len(primer) - primer.complementary_score <= 12,  # 18b
                                      primer.has_contig_gc_at(4, 5),  # 19
                                      primer.has_repeated_nucleotide(4),  # 20
                                      primer.gc < 0.40 or primer.gc > 0.60, # 21
                                      primer.contig_complementary_score >= 4,  # 22a
                                      len(primer) - primer.complementary_score <= 14,  # 22
                                      abs(len(primer)-22)))

def hamming(s1, s2):
    """Return the Hamming distance between equal-length sequences.
    Source: Wikipedia """
    s1 = str(s1)
    s2 = str(s2)
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length.")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))
