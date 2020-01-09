"""
License information goes here.
"""

from difflib import SequenceMatcher
import itertools

from Bio import Align
import regex

from .exceptions import StarpError

def add_tails(amas1, amas2, amplicon1, amplicon2, snp):
    """
    Args:
        amas1: An AmasPrimer object.
        amas2: An AmasPrimer object.
        amplicon1: The length of the amplicon from the first allele.
        amplicon2: The length of the amplicon from the second allele.
        snp: The snp these AmasPrimers were created around.
            This should be a substitution SNP.
    """
    from .models import Sequence

    tail1 = 'GCAACAGGAACCAGCTATGAC'
    tail2 = 'GACGCAAGTGAGCAGTATGAC'

    if amplicon1 - amplicon2 >= 8:
        amas1.tailed = Sequence(merge(tail1, str(amas1)))
        amas2.tailed = Sequence(merge(tail2, str(amas2)))
    elif amplicon1 - amplicon2 >= 1:
        amas1.tailed = Sequence(merge(tail2, str(amas1)))
        amas2.tailed = Sequence(merge(tail1, str(amas2)))
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

        if frozenset({str(amas1[0]), str(amas2[0])} == snp.nucleotides):
            if assigned_tail[snp.nucleotides][str(amas1[0])] == 1:
                # Assign tail1 to amas1, tail2 to amas2
                amas1.tailed = Sequence(merge(tail1, str(amas1)))
                amas2.tailed = Sequence(merge(tail2, str(amas2)))
            else:
                # Assign tail1 to amas2, tail2 to amas1
                amas1.tailed = Sequence(merge(tail2, str(amas1)))
                amas2.tailed = Sequence(merge(tail1, str(amas2)))
        elif frozenset({str(amas1[-1]), str(amas2[-1])} == snp.nucleotides):
            # Primers were created upstream
            if assigned_tail[snp.nucleotides][str(amas1[-1])] == 1:
                # Assign tail1 to amas1, tail2 to amas2
                amas1.tailed = Sequence(merge(tail1, str(amas1)))
                amas2.tailed = Sequence(merge(tail2, str(amas2)))
            else:
                # Assign tail1 to amas2, tail2 to amas1
                amas1.tailed = Sequence(merge(tail2, str(amas1)))
                amas2.tailed = Sequence(merge(tail1, str(amas2)))
        else:
            raise StarpError('Something went wrong when adding tails.')

    elif amplicon1 - amplicon2 >= -7:
        amas1.tailed = Sequence(merge(tail1, str(amas1)))
        amas2.tailed = Sequence(merge(tail2, str(amas2)))
    else:
        amas1.tailed = Sequence(merge(tail2, str(amas1)))
        amas2.tailed = Sequence(merge(tail1, str(amas2)))

    return (amas1, amas2)

def merge(str1: str, str2: str, max_chars=10000) -> str:
    """ Attempts to merge str1 and str2 up to max_chars characters.
    This is similar to concatentation, but an example works best.

    merge('ABCDEF', 'DEFGHI') -> 'ABCDEFGHI'

        ABCDEF
      +    DEFGHI
        ---------
        ABCDEFGHI

    """
    max_chars = min(len(str1), len(str2), max_chars)
    overlaps = 0

    for size in range(1, max_chars+1):
        print(str1[-1*size:], str2[:size])
        if str1[-1*size:] == str2[:size]:
            overlaps = size

    return str1 + str2[overlaps:]

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
    pattern = '(' + str(primer) + '){s<=' + str(TOLERANCE) + '}'

    matches_iter = iter([])
    for sequence in sequences:
        matches_iter = itertools.chain(matches_iter,
                                       regex.finditer(pattern,
                                                      str(sequence),
                                                      overlapped=True))

        for match in matches_iter:
            matched_seq = match.group()
            if (str(primer)[-1] == matched_seq[-1]
                    or hamming(str(primer)[-4:-1], matched_seq[-4:-1]) < 2):
                sites.append(match)

            if len(sites) >= stop:
                return sites

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

def record_spans(rprimers, allele1, allele2):
    """ Given a list of rprimers, record their span in both alleles.
    If it does not exist exactly in both alleles, remove it from the
    list. The span is recorded as the first match from regex .
    Any primers with multiple sites will be removed later. """
    updated_rprimers = list()

    for primer in rprimers:
        match = regex.search(str(primer.rev_comp()), str(allele1))
        print(str(primer.rev_comp()))
        if not match:
            continue

        primer.allele1_span = match.group(0)
        match = regex.search(str(primer.rev_comp()), str(allele2))

        if not match:
            continue

        primer.allele2_span = match.group(0)
        updated_rprimers.append(primer)

    return updated_rprimers

def rgenerate(ref_sequence, snp, min_length, max_length):
    """
    Return a list of all possible R primers from SNP end to
    the end of the reference sequence.

    Args:
        ref_sequence: The sequence to generate primers from.
        min_length: Minimum length of primer.
        max_length: Maximum length (inclusive) of primer.

    Returns:
        A list of all possible r primers.
    """
    from .models import Primer

    candidates = [Primer(str(ref_sequence[i:i+size].rev_comp()), i, i+size, -1)
                  for i in range(snp.position+1, len(ref_sequence))
                  for size in range(min_length, max_length+1)]

    return candidates

def rfilter_by_binding_sites(r_primers, allele1, allele2, nontargets,
                             max_num, amas):

    candidates = []

    for primer in r_primers:
        # The primer should have at most 1 binding site in each allele.
        # Primers are rev_comped to orient them the same as the alleles.
        if (len(binding_sites((allele1,), primer.rev_comp())) > 1
                or len(binding_sites((allele2,), primer.rev_comp())) > 1):
            continue

        # Move on if the primer has other potential binding sites in
        # these sequences.
        if binding_sites((allele1, allele2), primer, stop=1):
            continue

        if binding_sites(nontargets, primer, stop=1):
            continue

        if binding_sites(nontargets, primer.rev_comp(), stop=1):
            continue

        # Check complementarity with AMAS primers.
        if (len(primer) - complementary_score(primer.reverse(), amas[0]) > 5
                and len(primer) - complementary_score(primer.reverse(), amas[1]) > 5):
            candidates.append(primer)

        if len(candidates) >= max_num:
            break

    return candidates

def rfilter_by_binding_sites2(r_primers: list, sequences: tuple,
                              max_num: int, amas: tuple) -> list:
    """
    Filter R primers by checking their binding site on the sequences.
    One of these sequences should be the reference sequence, so one
    binding site should be guaranteed. However, extra binding sites on
    other sequences is strictly disallowed.

    If a primer has one binding site, the complementary scores between
    the primer and AMAS primers are calculated. If there are > 5
    mismatched nucleotides, add the primer to the list to be returned.
    At most, 3 primers will be returned.

    A binding site for a primer on the plus strand is defined to be a
    matching region with <= 4 nucleotide differences, with a 5'
    mismatch or < 2 mismatches at 2nd, 3rd, or 4th position from 5'
    end.

    A binding site for a primer on the minus strand is defined to be a
    matching region with <= 4 nucleotide differences, with a 3'
    mismatch or < 2 mismatches at 2nd, 3rd, or 4th position from 3'
    end.

    Both the primers and their reverse complements will be checked
    against the sequences.

    This is a very expensive operation, so the first 'max_num' good
    primers will be returned.

    Args:
        r_primers: The r_primers to check for binding sites.
        sequences: The sequences to check the primers against.
        max_num: Return the first 'max_num' primers with a single
            binding site. This helps with efficiency.
        amas: The 2-tuple containing the AMAS primers.

    Raises:
        StarpError: No binding sites could be found for a primer. All
            primers should have at least one binding site, since they
            were created from the reference sequence. This probably
            means the reference sequence is not one of the sequences.

    Returns:
        The primers with exactly one binding site.
    """

    one_binding_site_primers = list()

    for primer in r_primers:
        pattern = regex.compile('(' + str(primer) + '){s<=4}')
        rc_pattern = regex.compile('(' + str(primer.rev_comp()) + '){s<=4}')

        # Matches on the plus strand when the reverse primer is oriented the
        # same direction. This requires a reverse complement of the primer.
        rc_matches = iter([])

        # Matches on the sequences when the reverse primer is in its default
        # orientation.
        matches = iter([])

        for sequence in sequences:
            rc_matches = itertools.chain(rc_matches, regex.finditer(rc_pattern, str(sequence), overlapped=True))
            matches = itertools.chain(matches, regex.finditer(pattern, str(sequence), overlapped=True))

        # Remember, reverse primers are already on the minus strand. So, the
        # nucleotide differences on their 3' ends are checked.
        predicate = lambda match: (str(primer)[-1] == match.group()[-1]
                                   or hamming(str(primer)[-4:-1],
                                              match.group()[-4:-1]) < 2)
        binding_sites = filter(predicate, matches)


        rc_binding_sites = filter(lambda match: (str(primer.rev_comp())[0] == match.group()[0]
                                                 or hamming(str(primer.rev_comp())[1:4],
                                                            match.group()[1:4]) < 2),
                                  rc_matches)

        total_binding_sites = itertools.chain(binding_sites, rc_binding_sites)

        try:
            next(total_binding_sites)
        except StopIteration:
            # Weird... no binding sites could be found.
            raise StarpError('There was an issue with finding primer binding sites.')

        try:
            # If this succeeds, there are more than one binding sites.
            next(total_binding_sites)
        except StopIteration:
            # Primer has exactly one binding site.
            # Now check if it is complementary with the AMAS primers.
            if (len(primer) - complementary_score(primer.reverse(), amas[0]) > 5
                    and len(primer) - complementary_score(primer.reverse(), amas[1]) > 5):
                one_binding_site_primers.append(primer)

        if len(one_binding_site_primers) >= max_num:
            break

    return one_binding_site_primers

def rtailed(r_primers: list) -> list:
    """
    Add tails to reverse primers as described in
    docs/Starp R primer design[4311].

    Then, returns the primers shorter than 28 bases and those with
    tm <= 62 degrees C.

    Dr. Long adds {'', 'C', 'G', 'CG', 'GC', 'CGC', 'GCG'} to the ends
    of the primers, but he assumes they are oriented on the plus
    strand. However, reverse primers are already on the minus strand
    so these instead will be prepended to the primers.

    Args:
        r_primers: The primer list to add tails.

    Returns:
        A list of tailed R primers with length < 28 bp and tm <= 62
    """
    from .models import Primer

    tails = {'', 'C', 'G', 'CG', 'GC', 'CGC', 'GCG'}

    # The cartesian product of the tails and primers in a tuple.
    tails_and_primers = itertools.product(tails, r_primers)

    # Primers with tails.
    tailed = list()
    for tail, primer in tails_and_primers:
        if primer.start - len(tail) >= 0:
            tailed.append(Primer(tail + primer.sequence,
                                 primer.start-len(tail),
                                 primer.end, primer.strand))

    return [primer for primer in tailed
            if len(primer) < 28
            and primer.tm <= 62]


def rfilter(r_primers: list, amas: tuple, snps: list, pcr_max: int) -> list:
    """
    Removes the R primers meeting any of the following conditions:
    - overlap between primer and amas primer
    - Contains invalid characters. The only accepted alphabet
      is {A, C, G, T}.
    - The distance from amas[1].end and the primer is > pcr_max
    - A SNP exists in the primer region.
    - Has >= 10 contiguous G/C or >= 12 contiguous A/T
    - Has >= 8 As, Ts, Gs, or Cs
    - Has >= 6 di-nucleotide repeats
    - GC content is > 80% or < 20%

    Args:
        primers: Reverse primers, aka primers on -1 strand.
        amas: The 2-tuple of AMAS primers.
        snps: A list of SNP objects.
        pcr_max: An integer representing the longest allowed
            amplicon length.

    Returns:
        The primers that do not meet the above conditions.
    """

    candidates = [primer for primer in r_primers
                  if (not max(amas[0].end, amas[1].end) >= primer.start
                      and not regex.search('[^ACGT]', str(primer))
                      and abs(amas[1].end - primer.start) <= pcr_max
                      and not primer.has_contig_gc_at(10, 12)
                      and str(primer).count('A') < 8
                      and str(primer).count('T') < 8
                      and str(primer).count('G') < 8
                      and str(primer).count('C') < 8
                      and not primer.has_dinucleotide_repeat(6)
                      and 0.20 <= primer.gc <= 0.80)]

    candidates = rfilter_complementary(candidates)

    # Lastly, filter out the candidates with overlapping SNPs.
    for snp in snps:
        candidates = filter(lambda c: snp.position < c.start or snp.position >= c.end, candidates)

    return list(candidates)

def rfilter_complementary(r_primers: list) -> list:
    """
    Removes the R primers that have >= 10 contiguous self-complementary
    nucleotides or <= 4 nucleotides that are not complementary.
    """
    to_return = []
    for primer in r_primers:
        if (primer.contig_complementary_score < 10
                and len(primer) - primer.complementary_score > 5):
            to_return.append(primer)

    return to_return

    #return [primer for primer in r_primers
    #        if primer.contig_complementary_score < 10
    #        and len(primer) - primer.complementary_score > 5]

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

    Args:
        The R primers to be sorted.

    Returns:
        A sorted list of R primers.
    """

    # When a primer passes a condition, it should be evaluated to True
    # for that condition. This prevents a reversal of the list afterwards.
    return sorted(primers,
                  key=lambda primer: (not primer.has_contig_gc_at(9, 11),  # 1
                                      not primer.has_repeated_nucleotide(7),  # 2
                                      not primer.has_dinucleotide_repeat(5),  # 3
                                      0.25 <= primer.gc <= 0.75,  # 4
                                      primer.contig_complementary_score < 8,  # 5a
                                      len(primer) - primer.complementary_score > 6,  #5b
                                      not primer.has_contig_gc_at(8, 9),  # 6
                                      not primer.has_repeated_nucleotide(6),  # 7
                                      not primer.has_dinucleotide_repeat(4),  # 8
                                      0.30 <= primer.gc <= 0.70,  # 9
                                      primer.contig_complementary_score < 7,  # 10a
                                      len(primer) - primer.complementary_score > 8,  # 10b
                                      not primer.has_contig_gc_at(6, 7),  # 11
                                      not primer.has_repeated_nucleotide(5),  # 12
                                      not primer.has_dinucleotide_repeat(3),  # 13
                                      0.35 <= primer.gc <= 0.65,  # 14
                                      primer.contig_complementary_score < 6,  # 15a
                                      len(primer) - primer.complementary_score > 10,  # 15b
                                      not primer.has_in_last(5, 6, 7),  # 16
                                      not primer.has_in_last(3, 4, 4),  # 17
                                      primer.contig_complementary_score < 5,  # 18a
                                      len(primer) - primer.complementary_score > 12,  # 18b
                                      not primer.has_contig_gc_at(4, 5),  # 19
                                      not primer.has_repeated_nucleotide(4),  # 20
                                      0.40 <= primer.gc <= 0.60,  # 21
                                      primer.contig_complementary_score < 4,  # 22a
                                      len(primer) - primer.complementary_score > 14))  # 221

def hamming(s1, s2):
    """Return the Hamming distance between equal-length sequences.
    Source: Wikipedia """
    s1 = str(s1)
    s2 = str(s2)
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length.")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))
