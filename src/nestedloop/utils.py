from difflib import SequenceMatcher
import itertools

from Bio import Align
import regex

"""
This file provides a bunch of utility functions that are required in
other places.
"""

def binding_sites(sequences: tuple, primer, stop=2):
    """
    Returns a list of potential binding sites of this primer on these
    sequences. A binding site is defined as a match with <= 4 nucleotide
    subsitutions or a match on the 3' end or < 2 mismatches on
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
    """
    TOLERANCE = 4  # Number of acceptable nucleotide substitutions.

    sites = list()
    pattern = '(' + str(primer) + '){s<=' + str(TOLERANCE) + '}'
    revcomp_pattern = '(' + str(primer.rev_comp()) + '){s<=' + str(TOLERANCE) + '}'

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

            # If the final nucleotides are equal or there are 
            if (str(primer)[-1] == matched_seq[-1]
                    or hamming(str(primer)[-4:-1], matched_seq[-4:-1]) < 2):
                sites.append(match)

            if len(sites) >= stop:
                return sites

        matches_iter = regex.finditer(revcomp_pattern, str(sequence),
                                      overlapped=True)
        for match in matches_iter:
            if len(sites) >= stop:
                return sites

            matched_seq = match.group()
            if (str(primer)[-1] == matched_seq[-1]
                    or hamming(str(primer.rev_comp())[-4:-1], matched_seq[-4:-1]) < 2):
                sites.append(match)

        matches_iter = itertools.chain

    return sites

def complementary_score(s1: str, s2: str) -> int:
    """ Using Bio.Align from BioPython, calculate the score of a local
    sequence alignment using a very low open gap score as to completely
    eliminate gaps. Complement the reverse sequence because it is easier
    to check for same-ness than complementary nucleotides.

    Args:
        s1, s2: Strings of nucleotides.

    Returns:
        The maximum number of complementary nucleotides between s1 and
        s2 in a maximal alignment.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'

    # A low gap score makes sure that the sequence is not broken up into
    # smaller chunks in the alignment procedure.
    aligner.open_gap_score = -1000

    # score = the number of complementary nucleotides between s1 and s2.
    return aligner.score(str(s1), str(s2.complement()))

def contig_complementary_score(s1, s2) -> int:
    """ Complement the reverse sequence and check for equality. This is
    identical to checking for complementary of the reverse sequence.
    difflib.SequenceMatcher is very fast at this.

    Args:
        s1, s2: Sequence objects.

    Returns:
        The maximum number of contiguous complementary nucleotides
        between s1 and s2 in a maximal alignment.
    """

    matcher = SequenceMatcher()
    matcher.set_seq1(str(s1))
    matcher.set_seq2(str(s2.complement()))
    _, _, size = matcher.find_longest_match(0, len(s1), 0, len(s2))
    return size

def hamming(s1, s2):
    """
    Returns the hamming distance between s1 and s2. This is the
    number of mismatches between them.
    """
    assert len(s1) == len(s2)
    return sum(n1 != n2 for n1, n2 in zip(s1, s2))

def max_primer_length(tm_max):
    """ Returns the maximum primer length based on the maximum melting
    temperature. """
    return int((float(tm_max) + 7.5) / 2.5)

def min_primer_length(tm_min):
    """ Returns the minimum primer length based on the minimum melting
    temperature. """
    return int((float(tm_min) - 10) / 2.5)

def primer_generate(ref_sequence, strand: int, start: int, stop: int, tms: tuple):
    """
    Args:
        ref_sequence: The reference sequence to base primers from.
        strand: 1 corresponds to plus strand, -1 to minus strand.
        start: Primers are created in the range [start, stop].
        stop: See start.
        tms: A 3-tuple containing the (minimum, optimum, maximum)
            melting temperatures.

    Returns:
        A list of Primer candidates.

    Raises:
        ValueError: strand argument is not 1 or -1.
    """

    # This import at the top of the file causes circular imports.
    from .models import Primer

    if strand not in {1, -1}:
        raise ValueError(('Primer.strand must be a 1 or -1. Current value: '
                          '{primer.strand}.'))

    min_length = min_primer_length(tms[0])
    max_length = max_primer_length(tms[2])

    primers = []

    for i in range(start, stop):
        if strand == 1:
            # Generate all possible candidates with starting index equal
            # to i.
            candidates = [Primer(sequence=ref_sequence[i:i+size], span=(i, i+size), strand=strand)
                          for size in range(min_length, max_length+1)]
        else:
            # Generate all possible candidates with ending index equal to
            # i on the minus strand.
            candidates = [Primer(sequence=ref_sequence[i-size:i+1].rev_comp(), span=(i-size, i+1), strand=strand)
                          for size in range(min_length, max_length+1)]
            #candidates = [Primer(self.rev_comp_ref_sequence[i-size:i+1], i-size, i+1, strand=-1)
            #              for size in range(min_length, max_length+1)]

        # Choose the candidate with the optimum melting temperature.
        best_candidate = min(candidates, key=lambda c: abs(tms[1] - c.tm))

        # Verify that this primer is within the min/max temperatures and is
        # wholly contained within start/stop indices.
        if (tms[0] <= best_candidate.tm <= tms[2]
                and start <= best_candidate.start
                and stop >= best_candidate.end):
            primers.append(best_candidate)

    return primers

def pair_sort(pairs):
    """
    Returns a sorted list of pairs as defined in
    "docs/4_combine primers_NEW[2654].docx"

    The best pairs are at the lower indices.

    Args:
        pairs: The primer pairs to score.

    Returns:
        A sorted list of pairs.
    """

    return sorted(pairs,
                  key=lambda pair: (_complementary(pair, 4, 14),
                                    _complementary(pair, 5, 12),
                                    _complementary(pair, 6, 10),
                                    _complementary(pair, 7, 8),
                                    _complementary(pair, 8, 6),
                                    _complementary(pair, 9, 5)),
                  reverse=True)

def primer_sort(primers):
    """
    Returns a sorted list of primers as defined in
    "docs/nestedloop/3_design primers[2510].docx"

    The best primers are at lower indices.

    Note: Dr. Long's instructions say to sort reverse primers
    by 5' tacg. However, reverse primers have already been created
    5'->3' on the minus strand so this translates to 3' atgc
    because of the complementary nature of the minus strand.

    Args:
        primers: The list of primers to sort.

    Returns:
        A sorted list of primers.
    """

    return sorted(primers,
                  key=lambda primer: (primer.has_contig_gc_at(9,11),
                                      primer.has_mononucleotide_repeat(7),
                                      primer.has_dinucleotide_repeat(5),
                                      not 0.25 <= primer.gc <= 0.75,
                                      primer.has_contig_gc_at(8, 9),
                                      primer.has_mononucleotide_repeat(6),
                                      primer.has_dinucleotide_repeat(4),
                                      not 0.30 <= primer.gc <= 0.70,
                                      primer.has_contig_gc_at(6, 7),
                                      primer.has_mononucleotide_repeat(5),
                                      primer.has_dinucleotide_repeat(3),
                                      not 0.35 <= primer.gc <= 0.65,
                                      primer.has_in_last(5, 6, 7),
                                      primer.has_in_last(3, 4, 4),
                                      primer.has_contig_gc_at(4, 5),
                                      primer.has_mononucleotide_repeat(4),
                                      not 0.40 <= primer.gc <= 0.60,
                                      primer.is_self_complementary(8, 6),
                                      primer.is_self_complementary(7, 8),
                                      primer.is_self_complementary(6, 10),
                                      primer.is_self_complementary(5, 12),
                                      primer.is_self_complementary(4, 14),
                                      three_prime_atgc(primer)))

def three_prime_atgc(primer):
    """ Defines a mapping between {A, T, G, C} and integers that allows
    sorting with A > T > G > C. Here, we sort from the 3' end. """
    m = {'A': '0', 'T': '1', 'G': '2', 'C': '3'}
    seq = list(str(primer.reverse()))
    return ''.join([m[n] for n in seq])

def has_one_binding_site(primer, ref_sequence, nontarget_seqs):
    """ Searches for the primer sequence and its reverse complement
    through the reference sequence and all HSPs for any matches
    with <= 4 nucleotide differences. For each match, if the last
    nucleotide is a mismatch AND there are >= 2 mismatches in the
    2nd, 3rd, and 4th positions from the 3' end, then it is NOT
    considered a binding site. Else, increment num_binding_sites.

    Ex. The following IS a binding site.

        TGAGTCATGATCAGTATG
        |||  ||||| ||| |||
        TGAACCATGAACAGAATG

        But, this is NOT a binding site.

        TGAGTCATGATCAGTATG
        |||||||||| ||| |
        TGAGTCATGAACAGGACC

        Notice the mismatch on the 3' end. There are 3 mismatches
        in the last 4 bases.

    If more than one binding site is found, quit immediately.

    Args:
        primer: The primer for which we are searching fro binding
                sites.

    Returns:
        True if the primer has exactly one binding site. Otherwise,
        False.
    """
    sequences = [ref_sequence] + nontarget_seqs
    return len(binding_sites(sequences, primer)) == 1

###############
# Stuff that needs to be cleaned up.

def combine(f_primers, r_primers, num_to_return, pcr_min, pcr_max,
            ref_sequence, nontarget_seqs):
    """ Yunming's algorithm.
    f_primers and r_primers should be sorted such that the best primers
    are in position 0, 1, 2, etc. Then we step through the lists 20
    primers at a time and try combine them with primers at equal or lower
    index.

    Args:
        f_primers: A list of forward primers.
        r_primers: A list of reverse primers.
        num_to_return: The maximum number of pairs to return.

    Returns:
        The list of the best primer pairs.
    """
    from .models import Pair

    def update_used_primers(f, r):
        """ Add f and/or r to used_primers if custom primers were
        not set."""
        nonlocal used_primers
        if not f.custom:
            used_primers.update({f})
        if not r.custom:
            used_primers.update({r})

    # If both custom primers were given, then this is the only
    # combination we care about.
    if len(f_primers) == 1 and len(r_primers) == 1:
        if f_primers[0].custom and r_primers[0].custom:
            pairs = [Pair(f_primers[0], r_primers[0])]
            return pairs

    BLOCK_SIZE = 20
    used_primers = set()
    bad_primers = set()
    step = 0
    pairs = list()

    while (len(pairs) < num_to_return
               and (step < len(f_primers) or step < len(r_primers))):
        step += BLOCK_SIZE

        for primer in f_primers[step-BLOCK_SIZE:step] + r_primers[step-BLOCK_SIZE:step]:
            if not has_one_binding_site(primer, ref_sequence, nontarget_seqs):
                bad_primers.add(primer)

        # tuples is an iterator over tuples (f,r) with f in f_primers and r in r_primers.
        tuples = itertools.product(f_primers[:step],
                                   r_primers[step-BLOCK_SIZE:step])
        tuples = itertools.chain(tuples, itertools.product(f_primers[step-BLOCK_SIZE:step],
                                                           r_primers[:step]))

        for f, r in tuples:

            if (f in used_primers or f in bad_primers or r in used_primers or r in bad_primers):
                continue
            pair = Pair(f, r)
            if _good_combination(pair, pcr_min, pcr_max):
                pairs.append(pair)
                update_used_primers(f, r)

        pairs = pair_sort(pairs)

    return pairs[:num_to_return]

def _complementary(pair, score1, score2):
    """ Returns True if:
    1) The pair has contig complementary <= score1.
    2) Min Primer length - pair complementary >= score2.
    3) Each primer has contig comp score <= score1.
    """
    condition1 = pair.contig_complementary_score <= score1
    condition2 = (len(pair.forward_primer) - pair.complementary_score >= score2
                  and len(pair.reverse_primer) - pair.complementary_score >= score2)
    condition3 = (pair.forward_primer.contig_complementary_score <= score1
                  and pair.reverse_primer.contig_complementary_score <= score1)
    return condition1 and condition2 and condition3

def _good_combination(pair, pcr_min, pcr_max):
    """
    Return True if the following conditions are met:
    -The distance between the pairs is between pcr_min and pcr_max
    -The pair has <= 9 contig complementary
    -The length of each primer minus the pair complementary >= 5.
    """
    if not pcr_min <= pair.distance <= pcr_max:
        return False
    if pair.contig_complementary_score > 9:
        return False
    comp_score = pair.complementary_score
    if len(pair.forward_primer) - comp_score < 5:
        return False
    if len(pair.reverse_primer) - comp_score < 5:
        return False
    return True
