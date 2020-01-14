"""
License information goes here.
"""

from .parsers import TwoAlleles
from .exceptions import StarpError
from .models import Sequence, Snp, AmasPrimer

# G/C -> S and A/T -> W

# What the 4-tuple means:
# Let it be ('CG', -4, 'AT', -3)
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
        'SWPC' : -4, 'SWPG' : -3, 'SWNC' : ('CG', -4, 'AT', -3), 'SWNG' : ('CG', -4, 'AT', -3),
        'WSPC' : -4, 'WSPG' : -3, 'WSNC' : ('CG', -3, 'AT', -4), 'WSNG' : ('CG', -3, 'AT', -4),
        'WWPC' : -4, 'WWPG' : -3, 'WWNC' : -4, 'WWNG' : -3,
        'SPSC' : -4, 'SPSG' : -2, 'SPWC' : -4, 'SPWG' : -2,
        'SNSC' : -4, 'SNSG' : -2, 'SNWC' : ('CG', -4, 'AT', -2), 'SNWG' : ('CG', -4, 'AT', -2),
        'WPSC' : -4, 'WPSG' : -2, 'WPWC' : -4, 'WPWG' : -2,
        'WNSC' : ('CG', -2, 'AT', -4), 'WNSG' : ('CG', -2, 'AT', -4), 'WNWC' : -4, 'WNWG' : -2,
        'PSSC' : -3, 'PSSG' : -2, 'PSWC' : -3, 'PSWG' : -2,
        'PWSC' : -3, 'PWSG' : -2, 'PWWC' : -3, 'PWWG' : -2,
        'NSSC' : -3, 'NSSG' : -2, 'NSWC' : ('CG', -3, 'AT', -2), 'NSWG' : ('CG', -3, 'AT', -2),
        'NWSC' : ('CG', -2, 'AT', -3), 'NWSG' : ('CG', -2, 'AT', -3), 'NWWC' : -3, 'NWWG' : -2},
    frozenset(('C', 'T')) : {
        'SSPC' : -4, 'SSPT' : -3, 'SSNC' : -4, 'SSNT' : -3,
        'SWPC' : -4, 'SWPT' : -3, 'SWNC' : ('CG', -4, 'AT', -3), 'SWNT' : ('CG', -4, 'AT', -3),
        'WSPC' : -3, 'WSPT' : -4, 'WSNC' : ('CG', -3, 'AT', -4), 'WSNT' : ('CG', -3, 'AT', -4),
        'WWPC' : -4, 'WWPT' : -3, 'WWNC' : -4, 'WWNT' : -3,
        'SPSC' : -2, 'SPST' : -4, 'SPWC' : -4, 'SPWT' : -2,
        'SNSC' : -2, 'SNST' : -4, 'SNWC' : ('CG', -4, 'AT', -2), 'SNWT' : ('CG', -4, 'AT', -2),
        'WPSC' : -2, 'WPST' : -4, 'WPWC' : -2, 'WPWT' : -4,
        'WNSC' : ('CG', -2, 'AT', -2), 'WNST' : ('CG', -4, 'AT', -4), 'WNWC' : -2, 'WNWT' : -4,
        'PSSC' : -2, 'PSST' : -4, 'PSWC' : -3, 'PSWT' : -2,
        'PWSC' : -2, 'PWST' : -3, 'PWWC' : -2, 'PWWT' : -3,
        'NSSC' : -2, 'NSST' : -3, 'NSWC' : ('CG', -3, 'AT', -2), 'NSWT' : ('CG', -3, 'AT', -2),
        'NWGC' : ('CG', -2, 'AT', -2), 'NWGT' : ('CG', -3, 'AT', -3), 'NWWC' : -2, 'NWWT' : -3},
    frozenset(('C', 'A')) : {
        'SSPC' : -3, 'SSPA' : -4, 'SSNC' : -3, 'SSNA' : -4,
        'SWPC' : -4, 'SWPA' : -3, 'SWNC' : ('CG', -4, 'AT', -3), 'SWNA' : ('CG', -4, 'AT', -3),
        'WSPC' : -3, 'WSPA' : -4, 'WSNC' : ('CG', -3, 'AT', -3), 'WSNA' : ('CG', -4, 'AT', -4),
        'WWPC' : -3, 'WWPA' : -4, 'WWNC' : -3, 'WWNA' : -4,
        'SPSC' : -2, 'SPSA' : -4, 'SPWC' : -4, 'SPWA' : -2,
        'SNSC' : -2, 'SNSA' : -4, 'SNWC' : ('CG', -4, 'AT', -2), 'SNWA' : ('CG', -4, 'AT', -2),
        'WPSC' : -2, 'WPSA' : -4, 'WPWC' : -2, 'WPWA' : -4,
        'WNSC' : ('CG', -2, 'AT', -2), 'WNSA' : ('CG', -4, 'AT', -4), 'WNWC' : -2, 'WNWA' : -4,
        'PSSC' : -2, 'PSSA' : -3, 'PSWC' : -3, 'PSWA' : -2,
        'PWSC' : -2, 'PWSA' : -3, 'PWWC' : -2, 'PWWA' : -3,
        'NSSC' : -2, 'NSSA' : -3, 'NSWC' : ('CG', -3, 'AT', -2), 'NSWA' : ('CG', -3, 'AT', -2),
        'NWSC' : ('CG', -2, 'AT', -2), 'NWSA' : ('CG', -3, 'AT', -3), 'NWWC' : -2, 'NWWA' : -3},
    frozenset(('G', 'T')) : {
        'SSPT' : -3, 'SSPG' : -4, 'SSNT' : -3, 'SSNG' : -4,
        'SWPT' : -3, 'SWPG' : -4, 'SWNT' : ('CG', -4, 'AT', -3), 'SWNG' : ('CG', -4, 'AT', -3),
        'WSPT' : -4, 'WSPG' : -3, 'WSNT' : ('CG', -3, 'AT', -4), 'WSNG' : ('CG', -3, 'AT', -4),
        'WWPT' : -3, 'WWPG' : -4, 'WWNT' : -3, 'WWNG' : -4,
        'SPST' : -4, 'SPSG' : -2, '..W.' : -2, 'SPWG' : -4,
        'SNST' : -4, 'SNSG' : -2, 'SNWT' : ('CG', -4, 'AT', -2), 'SNWG' : ('CG', -4, 'AT', -2),
        'WPST' : -4, 'WPSG' : -2, 'WPWT' : -4, 'WPWG' : -2,
        'WNST' : ('CG', -4, 'AT', -4), 'WNSG' : ('CG', -2, 'AT', -2), 'WNWT' : -4, 'WNWG' : -2,
        'PSST' : -3, 'PSSG' : -2, 'PSWT' : -2, 'PSWG' : -3,
        'PWST' : -3, 'PWSG' : -2, 'PWWT' : -3, 'PWWG' : -2,
        'NSST' : -3, 'NSSG' : -2, 'NSWT' : ('CG', -3, 'AT', -2), 'NSWG' : ('CG', -3, 'AT', -2),
        'NWST' : ('CG', -3, 'AT', -3), 'NWSG' : ('CG', -2, 'AT', -2), 'NWWT' : -3, 'NWWG' : -2},
    frozenset(('G', 'A')) : {
        'SSPA' : -4, 'SSPG' : -3, 'SSNA' : -4, 'SSNG' : -3,
        'SWPA' : -3, 'SWPG' : -4, 'SWNA' : ('CG', -4, 'AT', -3), 'SWNG' : ('CG', -4, 'AT', -3),
        'WSPA' : -4, 'WSPG' : -3, 'WSNA' : ('CG', -3, 'AT', -4), 'WSNG' : ('CG', -3, 'AT', -4),
        'WWPA' : -4, 'WWPG' : -3, 'WWNA' : -4, 'WWNG' : -3,
        'SPSA' : -4, 'SPSG' : -2, 'SPWA' : -2, 'SPWG' : -4,
        'SNSA' : -4, 'SNSG' : -2, 'SNWA' : ('CG', -4, 'AT', -2), 'SNWG' : ('CG', -4, 'AT', -2),
        'WPSA' : -4, 'WPSG' : -2, 'WPWA' : -4, 'WPWG' : -2,
        'WNSA' : ('CG', -4, 'AT', -4), 'WNSG' : ('CG', -2, 'AT', -2), 'WNWA' : -4, 'WNWG' : -2,
        'PSSA' : -3, 'PSSG' : -2, 'PSWA' : -2, 'PSWG' : -3,
        'PWSA' : -3, 'PWSG' : -2, 'PWWA' : -3, 'PWWG' : -2,
        'NSSA' : -3, 'NSSG' : -2, 'NSWA' : ('CG', -3, 'AT', -2), 'NSWG' : ('CG', -3, 'AT', -2),
        'NWSA' : ('CG', -3, 'AT', -3), 'NWSG' : ('CG', -2, 'AT', -2), 'NWWA' : -3, 'NWWG' : -2},
    frozenset(('A', 'T')) : {
        'SSPA' : -4, 'SSPT' : -3, 'SSNA' : -4, 'SSNT' : -3,
        'SWPA' : -4, 'SWPT' : -3, 'SWNA' : ('CG', -4, 'AT', -3), 'SWNT' : ('CG', -4, 'AT', -3),
        'WSPA' : -4, 'WSPT' : -3, 'WSNA' : ('CG', -3, 'AT', -4), 'WSNT' : ('CG', -3, 'AT', -4),
        'WWPA' : -4, 'WWPT' : -3, 'WWNA' : -4, 'WWNT' : -3,
        'SPSA' : -4, 'SPST' : -2, 'SPWA' : -4, 'SPWT' : -2,
        'SNSA' : -4, 'SNST' : -2, 'SNWA' : ('CG', -4, 'AT', -2), 'SNWT' : ('CG', -4, 'AT', -2),
        'WPSA' : -4, 'WPST' : -2, 'WPWA' : -4, 'WPWT' : -2,
        'WNSA' : ('CG', -2, 'AT', -4), 'WNST' : ('CG', -2, 'AT', -4), 'WNWA' : -4, 'WNWT' : -2,
        'PSSA' : -3, 'PSST' : -2, 'PSWA' : -3, 'PSWT' : -2,
        'PWSA' : -3, 'PWST' : -2, 'PWWA' : -3, 'PWWT' : -2,
        'NSSA' : -3, 'NSST' : -2, 'NSWA' : ('CG', -3, 'AT', -2), 'NSWT' : ('CG', -3, 'AT', -2),
        'NWSA' : ('CG', -2, 'AT', -3), 'NWST' : ('CG', -2, 'AT', -3), 'NWWA' : -3, 'NWWT' : -2}}

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

def stretch_downstream(sequence, idx, length):
    """
    Dr. Long calls this technique "stretching", where you grab only the
    nucleotides and skip over any indel characters ('-'). For example,
    with sequence='A--GTACGGACT', idx=0, and length=4, we grab 'AGTA'.
    """

    important_seq = str(sequence[idx:]).replace('-', '')
    important_seq = important_seq.replace('-', '')
    try:
        stretched_seq = Sequence(important_seq[:length])
    except IndexError:
        raise StarpError('Unable to create Starp primers here.')

    return stretched_seq

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

def generate_amas_for_substitution(allele1, allele2, position, snps):
    """ Attempts to find the best upstream pair and the best
    downstream pair. These could be None.

    Does not substitute any bases.
    pos(ition) is relative to allele1, zero-indexed.
    """
    adj_position = convert_position(position, snps)

    pairs = zip(generate_amas_upstream(allele1, 1, position, 16, 26),
                generate_amas_upstream(allele2, 2, adj_position, 16, 26))
    upstream_pair = best_pair(pairs)

    pairs = zip(generate_amas_downstream(allele1, 1, position, 16, 26),
                generate_amas_downstream(allele2, 2, adj_position, 16, 26))
    downstream_pair = best_pair(pairs)

    if not upstream_pair and not downstream_pair:
        raise StarpError('Cannot find Starp primers at this location.')

    return upstream_pair, downstream_pair

def generate_amas_for_indel(allele1, allele2, position, snps):
    """ Does not substitute any bases.
    Based off of
    'how to design AMA-primers for Indel_20191125[3981].pptx' """
    adj_position = convert_position(position, snps)

    pairs = zip(generate_amas_upstream(allele1, 1, position, 1, 9),
                generate_amas_upstream(allele2, 2, adj_position, 1, 9))
    pairs = filter(lambda pair: pair[0][-1] != pair[1][-1], pairs)

    # Order the pairs based on
    # 1) The number of nucleotide differences in the last 4 bases.
    #    More differences is preferable.
    # 2) Length of the sequences in each pair. Longer is preferable.
    upstream_pair = sorted(
                        pairs,
                        key=lambda pair: (
                            Sequence.hamming(pair[0][-4:], pair[1][-4:]),
                            len(pair[0])),
                        reverse=True
                        )[0]

    pairs = zip(generate_amas_downstream(allele1, 1, position, 1, 9),
                generate_amas_downstream(allele2, 2, adj_position, 1, 9))

    # Remove the allele pairs having same base at 5' end
    pairs = filter(lambda pair: pair[0][0] != pair[1][0], pairs)

    downstream_pair = sorted(
                            pairs,
                            key=lambda pair: (
                                Sequence.hamming(pair[0][:4], pair[1][:4]),
                                len(pair[0])),
                            reverse=True
                            )[0]

    if not upstream_pair and not downstream_pair:
        raise StarpError('Cannot find Starp primers at this location.')

    return upstream_pair, downstream_pair

def generate_amas_upstream(allele, allele_num, idx, minimum, maximum):
    return [AmasPrimer(str(allele[idx-size:idx+1]), allele_num, (idx-size, idx+1), 'upstream')
            for size in range(minimum, maximum)
            if idx-size >= 0]

def generate_amas_downstream(allele, allele_num, idx, minimum, maximum):
    return [AmasPrimer(str(allele[idx:idx+size+1]), allele_num, (idx, idx+size+1), 'downstream')
            for size in range(minimum, maximum)
            if idx+size < len(allele)]

def substitute_bases(pair, snp, direction='upstream'):
    """
    Substitutes bases according to Dr. Long's instructions. This
    only works for substitutions and is undefined for insertions
    and deletions.

    The substitution principle he defines is:
    A -> C, T -> C, G -> A, C -> T

    I've tried my best to follow along to the instructions.

    To calculate the index that we need to substitute, all possible
    combinations are put into a dictionary. Since there are
    approximately 10*10*10*6*2 = 12,000 ways to arrange the last four
    nucleotides of an allele allowing multiple SNPs (the 10 being 4
    nucleotides plus 6 substitution SNPs, and the 2 being the number
    of nucleotides that define a SNP), there must be a way to
    reduce these combinations.

    The way this has been done is by grouping nucleotides and SNPs
    according to the following mapping:

    G->G
    C->G
    A->A
    T->A
    SNPs:
    (C/G) -> P  # Since C pairs with G and A to T, I'm mapping these
    (A/T) -> P  # to P for Paired.
    (C/A) -> N  # Since the next 4 SNPs contain nucleotides that do
    (C/T) -> N  # not pair with each other, I designate them N for
    (G/A) -> N  # Non-pair.
    (G/T) -> N

    Note that the last nucleotide is NOT mapped since its value
    must be known. For example, the last four nucleotides of an
    allele may be G(A/T)C(A/G). Using the mapping, this becomes
    GPG(A/G). Or, allele1_key = GPGA and allele2_key = GPGG. Then,
    the substitution index may be gleaned from the matries.

    This method works for both substitutions and indels.
    """

    if direction == 'downstream':
        # The SNP is at the beginning of the sequence, e.g.
        # pair[0] = XNNNNNNNN...
        # pair[1] = YNNNNNNNN...
        local_snps = TwoAlleles(f'>\n{pair[0][-4:-1]}\n>\n{pair[1][-4:-1]}').snps()
    elif direction == 'upstream':
        local_snps = TwoAlleles(f'>\n{pair[0][1:4]}\n>\n{pair[1][1:4]}').snps()

    if len(local_snps) == 0:
        new_amas1, new_amas2 = substitute_with_one_snp(pair, direction)
        pair[0].sequence = new_amas1
        pair[1].sequence = new_amas2
    elif len(local_snps) == 1:
        new_amas1, new_amas2 = substitute_with_two_snps(pair, snp, direction)
        pair[0].sequence = new_amas1
        pair[1].sequence = new_amas2

    return pair

def convert_position(pos, snps):
    """ Snp positions are given relative to allele1, so this function
    converts them into the corresponding position on allele2. """
    del_count = len([snp for snp in snps
                     if snp.position < pos and snp.type == 'insertion'])
    ins_count = len([snp for snp in snps
                     if snp.position < pos and snp.type == 'deletion'])

    return pos - del_count + ins_count

def substitute_with_one_snp(pair, direction='upstream'):
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

    if direction == 'downstream':
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

    if direction == 'downstream':
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    return (seq1, seq2)

def substitute_with_two_snps(pair, snp, direction='upstream'):
    """
    other_snp is the second snp towards the end that is not
    the final snp.

    Many times we only need to differentiate between a C/G or A/T snp,
    and all the rest. To do this, replace the Snp index with a P if it
    is a C/G or A/T snp (P for Paired) and N for all the rest.
    """
    pair = (str(pair[0]), str(pair[1]))

    if direction == 'downstream':
        pair = (pair[0][::-1], pair[1][::-1])

    if not pair[0][-1] != pair[1][-1]:
        raise ValueError('The sequences do not have a SNP in the last '
                         'position.')

    # SNPs at 2nd, 3rd, or 4th position from the 3' end.
    local_snps = TwoAlleles(f'>\n{pair[0][-4:-1]}\n>\n{pair[1][-4:-1]}').snps()

    if len(local_snps) != 1:
        raise ValueError('The sequences must have one SNP in the 2nd, 3rd, or '
                         '4th position from the 3\' end.')

    xsnp = local_snps[0]  # extra snp

    if xsnp.nucleotides == {'C', 'G'} or xsnp.nucleotides == {'A', 'T'}:
        placeholder = 'P'
    else:
        placeholder = 'N'

    seq1 = pair[0][:xsnp.position] + placeholder + pair[0][xsnp.position+1:]
    code = seq_to_ambiguity_code(pair[0][-4:-1]) + pair[0][-1]
    idx_to_sub = sub_index_two_snps[xsnp.nucleotides][code]
    seq1 = substitute(seq1, idx_to_sub)

    seq2 = pair[1][:xsnp.position] + placeholder + pair[1][xsnp.position+1:]
    code = seq_to_ambiguity_code(pair[1][-4:-1]) + pair[1][-1]
    idx_to_sub = sub_index_two_snps[xsnp.nucleotides][code]
    seq2 = substitute(seq2, idx_to_sub)

    if direction == 'downstream':
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

    return (seq1, seq2)
