"""
License information goes here.
"""

from .parsers import TwoAlleles
from .models import Sequence, Snp, AmasPrimer
from .utils import hamming

# These substitution matrices were created to deal with the
# combinatorial explosion of possible primer endings given by
# docs/starp/STARP F primer design[4312].docx
#
# Exempting the final nucleotide, 'G' and 'C' nucleotides were changed
# to an 'S' while 'A' and 'T' nucleotides were changed to a 'W'.
# 'P' represents a [C/G] OR [A/T] snps while 'N' represents all
# other substitution snps.
#
# A single integer means to substitute that index using the function
# substitute().
#
# A 2-tuple means to substitute those two indices.
#
# A 4-tuple is more complicated, and is only used when there is an
# additional SNP. An example describes its purpose best.
# Let it be ('S', -4, 'W', -3).
# Then, if the allele has a C or G in the additional SNP position,
# then substitute the -4th index. If it has an A or T in the
# additional SNP position, substitute the -3 index.

sub_index_one_snp = {
    frozenset(('C', 'G')) : {
        'SSSC' : -3, 'SSSG' : -2, 'SSWC' : (-3, -4), 'SSWG' : (-3, -4),
        'SWSC' : (-2, -4), 'SWSG' : (-2, -4), 'SWWC' : (-2, -3), 'SWWG' : (-2, -3),
        'WSSC' : (-2, -3), 'WSSG' : (-2, -3), 'WSWC' : (-2, -4), 'WSWG' : (-2, -4),
        'WWSC' : (-3, -4), 'WWSG' : (-3, -4), 'WWWC' : -3, 'WWWG' : -2},
    frozenset(('C', 'T')) : {
        'SSSC' : -2, 'SSST' : -3, 'SSWC' : -4, 'SSWT' : -3,
        'SWSC' : -2, 'SWST' : -4, 'SWWC' : -2, 'SWWT' : -3,
        'WSSC' : -2, 'WSST' : -3, 'WSWC' : -2, 'WSWT' : -4,
        'WWSC' : -4, 'WWST' : -3, 'WWWC' : -2, 'WWWT' : -3},
    frozenset(('C', 'A')) : {
        'SSSC' : -2, 'SSSA' : -3, 'SSWC' : -3, 'SSWA' : -4,
        'SWSC' : -2, 'SWSA' : -4, 'SWWC' : -2, 'SWWA' : -3,
        'WSSC' : -2, 'WSSA' : -3, 'WSWC' : -2, 'WSWA' : -4,
        'WWSC' : -3, 'WWSA' : -4, 'WWWC' : -2, 'WWWA' : -3},
    frozenset(('G', 'T')) : {
        'SSSG' : -2, 'SSST' : -3, 'SSWG' : -4, 'SSWT' : -3,
        'SWSG' : -2, 'SWST' : -4, 'SWWG' : -2, 'SWWT' : -4,
        'WSSG' : -2, 'WSST' : -3, 'WSWG' : -2, 'WSWT' : -4,
        'WWSG' : -4, 'WWST' : -3, 'WWWG' : -2, 'WWWT' : -2},
    frozenset(('G', 'A')) : {
        'SSSG' : -2, 'SSSA' : -3, 'SSWG' : -3, 'SSWA' : -4,
        'SWSG' : -2, 'SWSA' : -4, 'SWWG' : -2, 'SWWA' : -3,
        'WSSG' : -2, 'WSSA' : -3, 'WSWG' : -2, 'WSWA' : -4,
        'WWSG' : -3, 'WWSA' : -4, 'WWWG' : -2, 'WWWA' : -3},
    frozenset(('A', 'T')) : {
        'SSST' : -3, 'SSSA' : -4, 'SSWT' : -3, 'SSWA' : -4,
        'SWST' : -2, 'SWSA' : -4, 'SWWT' : -2, 'SWWA' : -3,
        'WSST' : -2, 'WSSA' : -3, 'WSWT' : -2, 'WSWA' : -4,
        'WWST' : -3, 'WWSA' : -4, 'WWWT' : -3, 'WWWA' : -4}}

sub_index_two_snps = {
    frozenset(('C', 'G')) : {
        'SSPC' : -4, 'SSPG' : -3, 'SSNC' : -4, 'SSNG' : -3,
        'SWPC' : -4, 'SWPG' : -3, 'SWNC' : ('S', -4, 'W', -3), 'SWNG' : ('S', -4, 'W', -3),
        'WSPC' : -4, 'WSPG' : -3, 'WSNC' : ('S', -3, 'W', -4), 'WSNG' : ('S', -3, 'W', -4),
        'WWPC' : -4, 'WWPG' : -3, 'WWNC' : -4, 'WWNG' : -3,
        'SPSC' : -4, 'SPSG' : -2, 'SPWC' : -4, 'SPWG' : -2,
        'SNSC' : -4, 'SNSG' : -2, 'SNWC' : ('S', -4, 'W', -2), 'SNWG' : ('S', -4, 'W', -2),
        'WPSC' : -4, 'WPSG' : -2, 'WPWC' : -4, 'WPWG' : -2,
        'WNSC' : ('S', -2, 'W', -4), 'WNSG' : ('S', -2, 'W', -4), 'WNWC' : -4, 'WNWG' : -2,
        'PSSC' : -3, 'PSSG' : -2, 'PSWC' : -3, 'PSWG' : -2,
        'PWSC' : -3, 'PWSG' : -2, 'PWWC' : -3, 'PWWG' : -2,
        'NSSC' : -3, 'NSSG' : -2, 'NSWC' : ('S', -3, 'W', -2), 'NSWG' : ('S', -3, 'W', -2),
        'NWSC' : ('S', -2, 'W', -3), 'NWSG' : ('S', -2, 'W', -3), 'NWWC' : -3, 'NWWG' : -2},
    frozenset(('C', 'T')) : {
        'SSPC' : -4, 'SSPT' : -3, 'SSNC' : -4, 'SSNT' : -3,
        'SWPC' : -4, 'SWPT' : -3, 'SWNC' : ('S', -4, 'W', -3), 'SWNT' : ('S', -4, 'W', -3),
        'WSPC' : -3, 'WSPT' : -4, 'WSNC' : ('S', -3, 'W', -4), 'WSNT' : ('S', -3, 'W', -4),
        'WWPC' : -4, 'WWPT' : -3, 'WWNC' : -4, 'WWNT' : -3,
        'SPSC' : -2, 'SPST' : -4, 'SPWC' : -4, 'SPWT' : -2,
        'SNSC' : -2, 'SNST' : -4, 'SNWC' : ('S', -4, 'W', -2), 'SNWT' : ('S', -4, 'W', -2),
        'WPSC' : -2, 'WPST' : -4, 'WPWC' : -2, 'WPWT' : -4,
        'WNSC' : ('S', -2, 'W', -2), 'WNST' : ('S', -4, 'W', -4), 'WNWC' : -2, 'WNWT' : -4,
        'PSSC' : -2, 'PSST' : -4, 'PSWC' : -3, 'PSWT' : -2,
        'PWSC' : -2, 'PWST' : -3, 'PWWC' : -2, 'PWWT' : -3,
        'NSSC' : -2, 'NSST' : -3, 'NSWC' : ('S', -3, 'W', -2), 'NSWT' : ('S', -3, 'W', -2),
        'NWSC' : ('S', -2, 'W', -2), 'NWST' : ('S', -3, 'W', -3), 'NWWC' : -2, 'NWWT' : -3},
    frozenset(('C', 'A')) : {
        'SSPC' : -3, 'SSPA' : -4, 'SSNC' : -3, 'SSNA' : -4,
        'SWPC' : -4, 'SWPA' : -3, 'SWNC' : ('S', -4, 'W', -3), 'SWNA' : ('S', -4, 'W', -3),
        'WSPC' : -3, 'WSPA' : -4, 'WSNC' : ('S', -3, 'W', -3), 'WSNA' : ('S', -4, 'W', -4),
        'WWPC' : -3, 'WWPA' : -4, 'WWNC' : -3, 'WWNA' : -4,
        'SPSC' : -2, 'SPSA' : -4, 'SPWC' : -4, 'SPWA' : -2,
        'SNSC' : -2, 'SNSA' : -4, 'SNWC' : ('S', -4, 'W', -2), 'SNWA' : ('S', -4, 'W', -2),
        'WPSC' : -2, 'WPSA' : -4, 'WPWC' : -2, 'WPWA' : -4,
        'WNSC' : ('S', -2, 'W', -2), 'WNSA' : ('S', -4, 'W', -4), 'WNWC' : -2, 'WNWA' : -4,
        'PSSC' : -2, 'PSSA' : -3, 'PSWC' : -3, 'PSWA' : -2,
        'PWSC' : -2, 'PWSA' : -3, 'PWWC' : -2, 'PWWA' : -3,
        'NSSC' : -2, 'NSSA' : -3, 'NSWC' : ('S', -3, 'W', -2), 'NSWA' : ('S', -3, 'W', -2),
        'NWSC' : ('S', -2, 'W', -2), 'NWSA' : ('S', -3, 'W', -3), 'NWWC' : -2, 'NWWA' : -3},
    frozenset(('G', 'T')) : {
        'SSPT' : -3, 'SSPG' : -4, 'SSNT' : -3, 'SSNG' : -4,
        'SWPT' : -3, 'SWPG' : -4, 'SWNT' : ('S', -4, 'W', -3), 'SWNG' : ('S', -4, 'W', -3),
        'WSPT' : -4, 'WSPG' : -3, 'WSNT' : ('S', -3, 'W', -4), 'WSNG' : ('S', -3, 'W', -4),
        'WWPT' : -3, 'WWPG' : -4, 'WWNT' : -3, 'WWNG' : -4,
        'SPST' : -4, 'SPSG' : -2, 'SPWT' : -2, 'SPWG' : -4,
        'SNST' : -4, 'SNSG' : -2, 'SNWT' : ('S', -4, 'W', -2), 'SNWG' : ('S', -4, 'W', -2),
        'WPST' : -4, 'WPSG' : -2, 'WPWT' : -4, 'WPWG' : -2,
        'WNST' : ('S', -4, 'W', -4), 'WNSG' : ('S', -2, 'W', -2), 'WNWT' : -4, 'WNWG' : -2,
        'PSST' : -3, 'PSSG' : -2, 'PSWT' : -2, 'PSWG' : -3,
        'PWST' : -3, 'PWSG' : -2, 'PWWT' : -3, 'PWWG' : -2,
        'NSST' : -3, 'NSSG' : -2, 'NSWT' : ('S', -3, 'W', -2), 'NSWG' : ('S', -3, 'W', -2),
        'NWST' : ('S', -3, 'W', -3), 'NWSG' : ('S', -2, 'W', -2), 'NWWT' : -3, 'NWWG' : -2},
    frozenset(('G', 'A')) : {
        'SSPA' : -4, 'SSPG' : -3, 'SSNA' : -4, 'SSNG' : -3,
        'SWPA' : -3, 'SWPG' : -4, 'SWNA' : ('S', -4, 'W', -3), 'SWNG' : ('S', -4, 'W', -3),
        'WSPA' : -4, 'WSPG' : -3, 'WSNA' : ('S', -3, 'W', -4), 'WSNG' : ('S', -3, 'W', -4),
        'WWPA' : -4, 'WWPG' : -3, 'WWNA' : -4, 'WWNG' : -3,
        'SPSA' : -4, 'SPSG' : -2, 'SPWA' : -2, 'SPWG' : -4,
        'SNSA' : -4, 'SNSG' : -2, 'SNWA' : ('S', -4, 'W', -2), 'SNWG' : ('S', -4, 'W', -2),
        'WPSA' : -4, 'WPSG' : -2, 'WPWA' : -4, 'WPWG' : -2,
        'WNSA' : ('S', -4, 'W', -4), 'WNSG' : ('S', -2, 'W', -2), 'WNWA' : -4, 'WNWG' : -2,
        'PSSA' : -3, 'PSSG' : -2, 'PSWA' : -2, 'PSWG' : -3,
        'PWSA' : -3, 'PWSG' : -2, 'PWWA' : -3, 'PWWG' : -2,
        'NSSA' : -3, 'NSSG' : -2, 'NSWA' : ('S', -3, 'W', -2), 'NSWG' : ('S', -3, 'W', -2),
        'NWSA' : ('S', -3, 'W', -3), 'NWSG' : ('S', -2, 'W', -2), 'NWWA' : -3, 'NWWG' : -2},
    frozenset(('A', 'T')) : {
        'SSPA' : -4, 'SSPT' : -3, 'SSNA' : -4, 'SSNT' : -3,
        'SWPA' : -4, 'SWPT' : -3, 'SWNA' : ('S', -4, 'W', -3), 'SWNT' : ('S', -4, 'W', -3),
        'WSPA' : -4, 'WSPT' : -3, 'WSNA' : ('S', -3, 'W', -4), 'WSNT' : ('S', -3, 'W', -4),
        'WWPA' : -4, 'WWPT' : -3, 'WWNA' : -4, 'WWNT' : -3,
        'SPSA' : -4, 'SPST' : -2, 'SPWA' : -4, 'SPWT' : -2,
        'SNSA' : -4, 'SNST' : -2, 'SNWA' : ('S', -4, 'W', -2), 'SNWT' : ('S', -4, 'W', -2),
        'WPSA' : -4, 'WPST' : -2, 'WPWA' : -4, 'WPWT' : -2,
        'WNSA' : ('S', -2, 'W', -4), 'WNST' : ('S', -2, 'W', -4), 'WNWA' : -4, 'WNWT' : -2,
        'PSSA' : -3, 'PSST' : -2, 'PSWA' : -3, 'PSWT' : -2,
        'PWSA' : -3, 'PWST' : -2, 'PWWA' : -3, 'PWWT' : -2,
        'NSSA' : -3, 'NSST' : -2, 'NSWA' : ('S', -3, 'W', -2), 'NSWT' : ('S', -3, 'W', -2),
        'NWSA' : ('S', -2, 'W', -3), 'NWST' : ('S', -2, 'W', -3), 'NWWA' : -3, 'NWWT' : -2}}


def differences_at_end(pair, snp_position):
    """
    If snp_position is 'last', return the number of nucleotide
    differences in the last 4 nucleotides. If snp_position is
    'first', return the number of differences in the first 4
    nucleotides.

    Args:
        pair: A 2-tuple of AMAS primers
        snp_position: Relative SNP position. Either 'first' or
            'last'.
    
    Returns:
        The hamming distance between the first or last 4 nucleotides.
    
    Raise:
        ValueError if either primer is shorter than 4 nucleotides.
    """

    seq1 = str(pair[0].sequence)
    seq2 = str(pair[1].sequence)

    if len(seq1) < 4 or len(seq2) < 4:
        raise ValueError('Primers must be greater than 4 nucleotides.')

    if snp_position == 'first':
        return hamming(seq1[:4], seq2[:4])
    else:
        return hamming(seq1[-4:], seq2[-4:])

def preserve_best_and_substitute(pairs, snp_position):
    """
    Corresponds to One SNP module, Two SNP module, and Three SNP module
    from STARP F primer design[4312].docx

    (1) Find the number of SNPs at the end of the shortest primer pair.
        This value influences the cutoff used when checking if the total
        SNP number between each pair is greater than 'cutoff'.

    (2) The temperature range and desired average temperature is chosen
        based on if there are any 'high snp pairs' from part 1. If there
        are none, then the substitute flag is set to true to signal that
        the selected primer pair will undergo substitutions.

    If no primers fit these conditions or the list of pairs is empty,
    None is returned.

    Args:
        pairs: A list of 2-tuples of AMAS pairs.
        snp_position: Position of the SNP relative to the AMAS primers.
            Should be either 'first' or 'last'.

    Returns:
        A 2-tuple representing the best AMAS pair.
        None if no primers have the described qualities.
    """

    # If the list is empty there is no good primer.
    if not pairs:
        return None

    # Flag for checking if substitutions need to be made on the selected
    # AMAS primer.
    substitute = False

    # Sort pairs by their length, shortest to longest.
    pairs = sorted(pairs, key=len)

    # Calculate SNP numbers at the ends of the first pair.
    snp_number = differences_at_end(pairs[0], snp_position)

    # All pairs should have 1 difference at the SNP location.
    if snp_number == 1:
        total_snp_number_cutoff = 4
    elif snp_number == 2:
        total_snp_number_cutoff = 3
    else:
        # All primers will be selected in this case.
        total_snp_number_cutoff = 0

    # Get pairs with >= 'total_snp_number_cutoff' SNPs
    high_snp_pairs = list(filter(
        lambda pair: hamming(pair[0].sequence, pair[1].sequence) >= total_snp_number_cutoff,
        pairs
    ))

    if high_snp_pairs:
        if snp_number > 2:
            low = 52
            high = 58
            average = 58
        else:
            # If there are some, try find the pair with Tm between 53
            # and 60 C and average is closest to 53 C.
            low = 53
            high = 60
            average = 53
    else:
        substitute = True
        low = 54
        high = 58
        average = 58

    # Keep pairs with Tm between 'low' and 'high'
    if high_snp_pairs:
        acceptable_pairs = list(filter(
            lambda pair: low <= pair[0].tm <= high and low <= pair[1].tm <= high,
            high_snp_pairs
        ))
    else:
        acceptable_pairs = list(filter(
            lambda pair: low <= pair[0].tm <= high and low <= pair[1].tm <= high,
            pairs
        ))

    # Sort pairs by average Tm closest to 'average' C.
    sorted_pairs = sorted(
        acceptable_pairs,
        key=lambda pair: abs(((pair[0].tm + pair[1].tm) / 2) - average)
    )

    # If there are no acceptable primer pairs, return None.
    if not sorted_pairs:
        return None

    # If there are acceptable primer pairs, select the first one from
    # the sorted list since it is closest to the desired temperature.
    best_pair = sorted_pairs[0]

    # If the flag is on, the bases of the best pair need to be
    # substituted.
    if substitute:
        best_pair[0].sequence, best_pair[1].sequence = substitute_bases((best_pair[0].sequence, best_pair[1].sequence), snp_position=snp_position)
    
    return best_pair

def amas_pair_filter(amas_pairs):
    """
    Removes primers pairs when one of the constituent primers has:
    1. >= 10 contiguous G/C or >= 12 contiguous A/T
    2. >= 8 of any single nucleotide
    3. >= 6 dinucleotide repeats
    4. GC > 0.80 or GC < 0.20

    Args:
        amas_pairs: A list of 2-tuples of AmasPrimers.
    """
    def filter_func(pair):
        for primer in pair:
            seq = Sequence(primer.sequence)
            if (seq.has_contig_gc_at(10, 12)
                    or seq.has_repeated_nucleotide(8)
                    or seq.has_dinucleotide_repeat(6)
                    or seq.gc < 0.20
                    or seq.gc > 0.80):
                return False
        return True

    return list(filter(filter_func, amas_pairs))

def substitute(seq, idx):
    """
    Substitutes the given index of the allele according to the
    mapping nucleotide_sub. This is necessary because strings
    and Sequence objects are immutable.

    Args:
        seq: the Sequence object to modify
        idx: an integer or tuple of the indices to change.

    Returns:
        The modified allele as a Sequence object.
    """
    if isinstance(idx, int):
        idx = (idx,)  # Change int to 1-tuple.

    seq = str(seq)

    nucleotide_sub = {'A' : 'C', 'T' : 'C', 'G' : 'A', 'C' : 'T'}
    seq = list(seq)

    for i in idx:
        seq[i] = nucleotide_sub[seq[i]]

    return Sequence(''.join(seq))

def seq_to_ambiguity_code(sequence: str):
    """ Converts 'C' and 'G' to 'S', and 'A'/'T' to 'W' as defined by
    http://www.reverse-complement.com/ambiguity.html


    Note that these are not being used to designate SNPs. In Dr. Long's
    instructions, many times he asks for a count of how many G/C bases
    are in a certain sequence, and this function reduces the sizes of
    the already massive dictionaries above.
    """
    return (sequence.replace('C', 'S').replace('G', 'S')
            .replace('A', 'W').replace('T', 'W'))

def generate_amas_for_substitution(allele1, allele2, position):
    """ Attempts to find the best upstream pair and the best
    downstream pair. These could be None.
    
    Args:
        allele1: The aligned first allele.
        allele2: The aligned second allele.
        position: The index around which to generate primers.
    """
    pairs = list(zip(generate_amas_upstream(allele1, 1, position, 16, 26),
                     generate_amas_upstream(allele2, 2, position, 16, 26)))

    # Remove pairs whose primers have undesirable characteristics.
    pairs = amas_pair_filter(pairs)

    upstream_pair = preserve_best_and_substitute(pairs, snp_position='last')

    pairs = list(zip(generate_amas_downstream(allele1, 1, position, 16, 26),
                     generate_amas_downstream(allele2, 2, position, 16, 26)))

    pairs = amas_pair_filter(pairs)
    
    downstream_pair = preserve_best_and_substitute(pairs, snp_position='first')

    return upstream_pair, downstream_pair

def generate_amas_for_indel(allele1, allele2, position):
    """
    
    An example works best.

    Say we have the alleles

        GTGG ACGCTCGAGGACTATAG--TCAGGAGAGGTGGGCATGG
        |||| |||||||||||||||||  |||||||||||||||||||
        GTGG ACGCTCGAGGACTATAGTCTCAGGAGAGGTGGGCATGG

    Then the upstream pairs that get generated are

        ACGCTCGAGGACTATAGT
        ||||||||||||||||||
        ACGCTCGAGGACTATAGT

        ACGCTCGAGGACTATAGTC
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTC

        ACGCTCGAGGACTATAGTCA
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCT

        ACGCTCGAGGACTATAGTCAG
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCTC

        ACGCTCGAGGACTATAGTCAGG
        |||||||||||||||||||
        ACGCTCGAGGACTATAGTCTCA

        ...

        ACGCTCGAGGACTATAGTCAGGAGA
        |||||||||||||||||||    ||
        ACGCTCGAGGACTATAGTCTCAGGA

    Args:
        allele1: The aligned first allele.
        allele2: The aligned second allele.
        position: The index of the indel around which to generate
            primers.

    Returns:
        upstream_pair, downstream_pair

        upstream_pair: The best AMAS pair where the majority of the
            primer sequence is upstream of the SNP.
        downstream_pair: The best AMAS pair where the majority of the
            primer sequence is downstream from the SNP.
    """

    # Create upstream AMAS primers by moving the SNP back 15 positions
    # and creating downstream primers from that new position.
    pairs = zip(generate_amas_downstream(allele1, 1, position-15, 16, 26),
                     generate_amas_downstream(allele2, 2, position-15, 16, 26))

    # Remove pairs with the same nucleotide on the 3' end.
    pairs = list(filter(lambda pair: pair[0][-1] != pair[1][-1], pairs))

    # Filter pairs for undesirable characteristics.
    pairs = amas_pair_filter(pairs)

    # Order the pairs based on
    # 1) The number of nucleotide differences in the last 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Shorter is preferable.
    upstream_pairs = sorted(pairs,
                            key=lambda pair: (
                                Sequence.hamming(pair[0][-4:], pair[1][-4:]),
                                1 / len(pair[0])),
                            reverse=True)

    # Go through the list, verifying the melting temperatures are good and
    # substituting bases if needed. Since these are already sorted, the first
    # primer that succeeds is the one sought after.
    upstream_pair = None
    for pair in upstream_pairs:
        upstream_pair = preserve_best_and_substitute([pair], snp_position='last')
        if upstream_pair:
            break

    # Create downstream AMAS primers by moving the SNP forward 15
    # positions and creating upstream primers from that position.
    pairs = zip(generate_amas_upstream(allele1, 1, position+15, 16, 26),
                generate_amas_upstream(allele2, 2, position+15, 16, 26))

    # Remove pairs with the same nucleotide on the 5' end.
    pairs = list(filter(lambda pair: pair[0][0] != pair[1][0], pairs))

    # Filter pairs for undesirable characteristics
    pairs = amas_pair_filter(pairs)

    # Order the pairs based on
    # 1) The number of nucleotide differences in the first 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Shorter is preferable.
    downstream_pairs = sorted(pairs,
                              key=lambda pair: (
                                  Sequence.hamming(pair[0][:4], pair[1][:4]),
                                  1 / len(pair[0])),
                              reverse=True)

    # Go through the list, verifying the melting temperatures are good and
    # substituting bases if needed. Since these are already sorted, the first
    # primer that succeeds is the one sought after.
    downstream_pair = None
    for pair in downstream_pairs:
        downstream_pair = preserve_best_and_substitute([pair], snp_position='first')
        if downstream_pair:
            break

    return upstream_pair, downstream_pair

def generate_amas_upstream(allele, allele_num, pos, minimum, maximum):
    """ Returns AMAS primers upstream of the position using the 
    aligned sequence. 

    Args:
        allele: The aligned allele to create primers from. It is
            possible for this to contain dashes.
        allele_num: This is using either allele 1 or 2. This is
            necessary for the primer instantiation.
        pos: The position to start from.
        minimum: The minimum primer length.
        maximum: The maximum primer length, inclusive.
    """
    sliced = str(allele[:pos+1]).replace('-', '')
    return [AmasPrimer(sliced[0-size:],
                       allele_num,
                       span=(len(sliced)-size, len(sliced)))
            for size in range(minimum, maximum+1)
            if pos >= size]

def generate_amas_downstream(allele, allele_num, pos, minimum, maximum):
    """ Returns a list of AMAS primers downstream of the position using
    the aligned allele. Doing this downstream is slightly more
    complicated since we have to consider the positions of the dashes.
    """
    dashes_before_pos = str(allele[:pos]).count('-')
    dashes_after_pos = str(allele[pos+1:]).count('-')
    is_dash = 1 if allele[pos] == '-' else 0

    sliced = str(allele[pos:]).replace('-', '')
    return [AmasPrimer(sliced[:size],
                       allele_num,
                       span=(pos-dashes_before_pos,
                             pos-dashes_before_pos+size))
            for size in range(minimum, maximum+1)
            if size <= len(sliced)]

def substitute_bases(pair, snp_position='last'):
    """
    Substitutes bases according to Dr. Long's instructions.
    See docs/STARP F primer design[4312].docx, pages 2-14.

    Args:
        pair: A 2-tuple of Sequence objects or strings.
        snp_position: should be either 'first' or 'last', signifying
            if the snp is at the beginning or end of the pair.
    """
    if not pair:
        return None

    if not ((type(pair[0]) == Sequence and type(pair[1]) == Sequence)
            or (type(pair[0]) == str and type(pair[1]) == str)):
        raise ValueError('Pair must be a tuple of strings or Sequences')

    # All of the keys in the large dicts require 4+ characters.
    if len(pair[0]) < 4:
        return pair

    if snp_position == 'first':
        str_pair = (str(Sequence(pair[0]).rev_comp()), str(Sequence(pair[1]).rev_comp()))
    else:
        str_pair = (str(pair[0]), str(pair[1]))

    snp = Snp(f'.{len(str_pair[0])-1}{str_pair[0][-1]}>{str_pair[1][-1]}')
    local_snps = TwoAlleles(f'>\n{str_pair[0][-4:-1]}\n>\n{str_pair[1][-4:-1]}').snps()

    if len(local_snps) == 0:
        new_amas1, new_amas2 = substitute_with_one_snp(str_pair, 'last')
        pair = (str(new_amas1), str(new_amas2))
    elif len(local_snps) == 1:
        new_amas1, new_amas2 = substitute_with_two_snps(str_pair, snp, 'last')
        pair = (str(new_amas1), str(new_amas2))

    pair = (Sequence(pair[0]), Sequence(pair[1]))

    # Reorient primers
    if snp_position == 'first':
        pair = (pair[0].rev_comp(), pair[1].rev_comp())

    return pair

def substitute_with_one_snp(pair, snp_position='last'):
    """
    Substitute the bases of the pair's sequences when the only SNP
    occurs at either the first or last nucleotide. This function probably
    should only be called by substitute_bases().

    This function assumes the primers in the pair are oriented on the
    plus strand.

    Notation for ambiguity codes comes from
    http://www.reverse-complement.com/ambiguity.html

    Ambiguity codes:
    G/C = S
    A/T = W

    Args:
        pair: A tuple of same-length strings or sequences.
        snp_position: Acceptable values are 'first' and 'last'.

    Returns:
        A 2-tuple of the pair with appropriate bases substituted
            of type Sequence.
    """

    # If the pair are downstream primers, as noted by snp_position
    # being 'first', we can reverse complement them and still use
    # the same instructions.
    if snp_position == 'first':
        pair = (str(pair[0].rev_comp()), str(pair[1].rev_comp()))
    else:
        pair = (str(pair[0]), str(pair[1]))

    if not pair[0][-1] != pair[1][-1]:
        raise ValueError('The sequences do not have a SNP in the last '
                         'position.')

    if pair[0][-4:-1] != pair[1][-4:-1]:
        raise ValueError('The sequences must be equal in the 2nd, 3rd, and '
                         '4th position from the 3\' end.')

    snp = Snp(f'.{len(pair[0])}{pair[0][-1]}>{pair[1][-1]}')

    # The only Snp between the two sequences is at the last index.
    code = seq_to_ambiguity_code(pair[0][-4:-1]) + pair[0][-1]
    idx_to_sub = sub_index_one_snp[snp.nucleotides][code]
    seq1 = substitute(pair[0], idx_to_sub)

    code = seq_to_ambiguity_code(pair[1][-4:-1]) + pair[1][-1]
    idx_to_sub = sub_index_one_snp[snp.nucleotides][code]
    seq2 = substitute(pair[1], idx_to_sub)

    # Orient the sequences back to their original orientation.
    if snp_position == 'first':
        seq1 = seq1.rev_comp()
        seq2 = seq2.rev_comp()

    return (seq1, seq2)

def substitute_with_two_snps(pair, snp, snp_position='last'):
    """
    Substitute the bases of the pair sequences when there are two SNPs
    between the sequences in the last four bases. This function probably
    should only be called by substitute_bases().

    This function assumes the primers in the pair are oriented on the
    plus strand.

    For reference, see
    docs/STARP F primer design[4312].docx
    """
    # Long's written instructions assume the SNP is placed at the end of
    # the sequences. However, this function will also be called with a
    # SNP at the beginning of the sequences. In this case, his
    # verbal instructions were to perform the same operations but at the
    # beginning of the sequences. To avoid making another very similar
    # function, the sequences are reversed here then returned to their
    # original orientation at the end of the function.

    # If the pair are downstream primers, as noted by snp_position
    # being 'first', we can reverse complement them and still use
    # the same instructions.
    if snp_position == 'first':
        pair = (str(pair[0].rev_comp()), str(pair[1].rev_comp()))
    else:
        pair = (str(pair[0]), str(pair[1]))

    # Long's written instructions assume the SNP is placed at the end of
    # the sequences. However, this function will also be called with a
    # SNP at the beginning of the sequences. In this case, his
    # verbal instructions were to perform the same operations but at the
    # beginning of the sequences. To avoid making another very similar
    # function, the sequences are reversed here then returned to their
    # original orientation at the end of the function.

    if not pair[0][-1] != pair[1][-1]:
        raise ValueError('The sequences do not have a SNP in the last '
                         'position.')

    # All SNPs between the two alleles.
    snps = TwoAlleles(f'>\n{pair[0]}\n>\n{pair[1]}').snps()

    # SNPs at 2nd, 3rd, or 4th position from the 3' end.
    snps = list(filter(lambda snp: 0 < len(pair[0])-snp.position-1 < 4, snps))

    if len(snps) != 1:
        raise ValueError('The sequences must have one SNP in the 2nd, 3rd, or '
                         '4th position from the 3\' end.')

    xsnp = snps[0]  # extra snp

    # Part of the codes created later concerns themselves about whether
    # the additional SNP is a [C/G] or [A/T] SNP, or the others (see
    # STARP F primer design[4312].docx, page 6). To simplify the possible
    # number of codes, [C/G] and [A/T] SNPs are designated as 'P' for paired
    # and all the rest are designated 'N'.
    if xsnp.nucleotides == {'C', 'G'} or xsnp.nucleotides == {'A', 'T'}:
        placeholder = 'P'
    else:
        placeholder = 'N'

    encoded_seq = pair[0][:xsnp.position] + placeholder + pair[0][xsnp.position+1:]
    code = seq_to_ambiguity_code(encoded_seq[-4:-1]) + pair[0][-1]
    idx_to_sub = sub_index_two_snps[snp.nucleotides][code]

    # Check if a 4-tuple was returned from the dict.
    if isinstance(idx_to_sub, tuple):
        if len(idx_to_sub) == 4:
            if encoded_seq[xsnp.position] == idx_to_sub[0]:
                idx_to_sub = idx_to_sub[1]
            else:
                idx_to_sub = idx_to_sub[3]

    seq1 = substitute(pair[0], idx_to_sub)

    encoded_seq = pair[1][:xsnp.position] + placeholder + pair[1][xsnp.position+1:]
    code = seq_to_ambiguity_code(encoded_seq[-4:-1]) + pair[1][-1]
    idx_to_sub = sub_index_two_snps[snp.nucleotides][code]

    # Check if a 4-tuple was returned from the dict.
    if isinstance(idx_to_sub, tuple):
        if len(idx_to_sub) == 4:
            if encoded_seq[xsnp.position] == idx_to_sub[0]:
                idx_to_sub = idx_to_sub[1]
            else:
                idx_to_sub = idx_to_sub[3]

    seq2 = substitute(pair[1], idx_to_sub)

    # Orient the sequences back to their original orientation.
    if snp_position == 'first':
        seq1 = seq1.rev_comp()
        seq2 = seq2.rev_comp()

    return (seq1, seq2)
