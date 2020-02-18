"""
License information goes here.
"""

from .parsers import TwoAlleles
from .exceptions import StarpError
from .models import Sequence, Snp, AmasPrimer

# G/C -> S and A/T -> W

# What the 4-tuple means:
# Let it be ('S', -4, 'W', -3)
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
        'NWGC' : ('S', -2, 'W', -2), 'NWGT' : ('S', -3, 'W', -3), 'NWWC' : -2, 'NWWT' : -3},
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

def best_pair(pairs):
    """
    Chooses the best AMAS primer pair.
    The best pair is the one with an average melting temperature
    closest to 58 C, with both primer's Tm between 54 and 58 C.

    If no pairs have this quality, then return None.

    Args:
        pairs: A list of 2-tuples of AMAS pairs.

    Returns:
        A 2-tuple representing the best AMAS pair if a pair exists
        with both primers having a Tm between 54 and 58 C.

        None if a pair with this quality does not exist.

    Raises:
        None.
    """
    tm_dict = dict()

    # Put pairs and melting temperatures in a dict.
    # (pair[0], pair[1]) : (tm1, tm2)
    tm_dict = {pair : (pair[0].tm, pair[1].tm) for pair in pairs}

    # Keep the pairs with both melting temperatures between 54 and 58 C.
    good_pairs = {pair : tm_dict[pair]
                  for pair in tm_dict
                  if (54 <= tm_dict[pair][0] <= 58
                      and 54 <= tm_dict[pair][1] <= 58)}

    # Compute the average melting temperatures between
    average_tms = {pair : (good_pairs[pair][0]+good_pairs[pair][1])/2
                   for pair in good_pairs}

    best_pairs = sorted(average_tms, key=lambda x: abs(average_tms[x]-58))

    if best_pairs == []:
        best_pair = None
    else:
        best_pair = best_pairs[0]

    return best_pair

def substitute(allele, idx):
    """
    Substitutes the given index of the allele according to the
    mapping nucleotide_sub. This is necessary because strings
    and Sequence objects are immutable.

    Args:
        allele: the Sequence object to modify
        idx: an integer or tuple of the indices to change.

    Returns:
        The modified allele as a Sequence object.
    """
    if isinstance(idx, int):
        idx = (idx,)  # Change int to 1-tuple.

    allele = str(allele)

    nucleotide_sub = {'A' : 'C', 'T' : 'C', 'G' : 'A', 'C' : 'T'}
    allele = list(allele)

    for i in idx:
        allele[i] = nucleotide_sub[allele[i]]

    return Sequence(''.join(allele))


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

    Does not substitute any bases.
    
    Args:
        allele1: The aligned first allele.
        allele2: The aligned second allele.
        position: The index around which to generate primers.
    """
    pairs = zip(generate_amas_upstream(allele1, 1, position, 16, 26),
                generate_amas_upstream(allele2, 2, position, 16, 26))
    upstream_pair = best_pair(pairs)

    pairs = zip(generate_amas_downstream(allele1, 1, position, 16, 26),
                generate_amas_downstream(allele2, 2, position, 16, 26))
    downstream_pair = best_pair(pairs)

    if not upstream_pair and not downstream_pair:
        raise StarpError('Cannot find Starp primers at this location.')

    return upstream_pair, downstream_pair

def generate_amas_for_indel(allele1, allele2, position):
    """ Does not substitute any bases.
    
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

    # Create upstream AMAS primers by moving the SNP back 17 positions
    # and creating downstream primers from that new position.
    pairs = zip(generate_amas_downstream(allele1, 1, position-17, 18, 26),
                generate_amas_downstream(allele2, 2, position-17, 18, 26))

    # Remove pairs with the same nucleotide on the 3' end.
    pairs = filter(lambda pair: pair[0][-1] != pair[1][-1], pairs)

    # Order the pairs based on
    # 1) The number of nucleotide differences in the last 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Longer is preferable.
    upstream_pairs = sorted(pairs,
                            key=lambda pair: (
                                Sequence.hamming(pair[0][-4:], pair[1][-4:]),
                                len(pair[0])),
                            reverse=True
                            )

    upstream_pair = upstream_pairs[0] if upstream_pairs else []

    # Create downstream AMAS primers by moving the SNP forward 17
    # positions and creating upstream primers from that position.
    pairs = zip(generate_amas_upstream(allele1, 1, position+17, 18, 26),
                generate_amas_upstream(allele2, 2, position+17, 18, 26))

    # Remove pairs with the same nucleotide on the 5' end.
    pairs = filter(lambda pair: pair[0][0] != pair[1][0], pairs)

    # Order the pairs based on
    # 1) The number of nucleotide differences in the first 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Longer is preferable.
    downstream_pairs = sorted(
                             pairs,
                             key=lambda pair: (
                                 Sequence.hamming(pair[0][:4], pair[1][:4]),
                                 len(pair[0])),
                             reverse=True
                             )

    downstream_pair = downstream_pairs[0] if downstream_pairs else []

    if not upstream_pair and not downstream_pair:
        raise StarpError('Cannot find Starp primers at this location.')

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

    snp_position should be either 'first' or 'last', signifying
    if the snp is at the beginning or end of the pair.
    """
    if not pair:
        return None

    # All of the keys in the large dicts require 4+ characters.
    if len(pair[0]) < 4:
        return pair

    if snp_position == 'first':
        # The SNP is at the beginning of the sequence, e.g.
        # pair[0] = XNNNNNNNN...
        # pair[1] = YNNNNNNNN...
        snp = TwoAlleles(f'>\n{pair[0][0]}\n>\n{pair[1][0]}').snps()[0]
        local_snps = TwoAlleles(f'>\n{pair[0][1:4]}\n>\n{pair[1][1:4]}').snps()
    elif snp_position == 'last':
        snp = TwoAlleles(f'>\n{pair[0][-1]}\n>\n{pair[1][-1]}').snps()[0]
        local_snps = TwoAlleles(f'>\n{pair[0][-4:-1]}\n>\n{pair[1][-4:-1]}').snps()

    if len(local_snps) == 0:
        new_amas1, new_amas2 = substitute_with_one_snp(pair, snp_position)
        pair[0].sequence = new_amas1
        pair[1].sequence = new_amas2
    elif len(local_snps) == 1:
        new_amas1, new_amas2 = substitute_with_two_snps(pair, snp, snp_position)
        pair[0].sequence = new_amas1
        pair[1].sequence = new_amas2

    return pair



def substitute_with_one_snp(pair, snp_position='last'):
    """
    http://www.reverse-complement.com/ambiguity.html
    Ambiguity codes:
    G/C = S
    A/T = W

    Args:
        pair: A tuple of same-length strings or sequences.

    Returns:
        A 2-tuple of the pair with appropriate bases substituted
            of type Sequence.
    """

    pair = (str(pair[0]), str(pair[1]))

    if snp_position == 'first':
        pair = (pair[0][::-1], pair[1][::-1])

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

    if snp_position == 'first':
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    return (seq1, seq2)

def substitute_with_two_snps(pair, snp, snp_position='last'):
    """
    Substitute the bases of the pair sequences when there are two SNPs
    between the sequences in the last four bases.

    For reference, see
    docs/STARP F primer design[4312].docx
    """
    pair = (str(pair[0]), str(pair[1]))

    # Long's written instructions assume the SNP is placed at the end of
    # the sequences. However, this function will also be called with a
    # SNP at the beginning of the sequences. In this case, his
    # verbal instructions were to perform the same operations but at the
    # beginning of the sequences. To avoid making another very similar
    # function, the sequences are reversed here then returned to their
    # original orientation at the end of the function.
    if snp_position == 'first':
        pair = (pair[0][::-1], pair[1][::-1])

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
    if len(idx_to_sub) == 4:
        if encoded_seq[xsnp.position] == idx_to_sub[0]:
            idx_to_sub = idx_to_sub[1]
        else:
            idx_to_sub = idx_to_sub[3]

    seq2 = substitute(pair[1], idx_to_sub)

    # Return the sequences to their original orientation if they were
    # reversed at the start of the method.
    if snp_position == 'first':
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    return (seq1, seq2)
