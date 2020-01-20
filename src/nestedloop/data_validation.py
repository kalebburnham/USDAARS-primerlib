import regex

from .models import Primer, Sequence
from .parsers import PairwiseParser, XmlParser
from .utils import binding_sites

def validate_reference_sequence(ref_sequence: str):
    """ Removes the FASTA header (if there is one) from the input which is
    assumed to be on the first line then strips any new line characters or
    whitespace at the beginning and end. """
    from nestedloop import REF_SEQ_MIN_LENGTH, REF_SEQ_MAX_LENGTH

    if not isinstance(ref_sequence, str):
        raise TypeError('Reference sequence must be a string.')

    if ref_sequence.startswith('>'):
        ref_sequence = ''.join(ref_sequence.splitlines()[1:])

    ref_sequence = ref_sequence.strip().upper()

    if not REF_SEQ_MIN_LENGTH <= len(ref_sequence) <= REF_SEQ_MAX_LENGTH:
        raise ValueError(('Reference sequence must contain between '
                          + str(REF_SEQ_MIN_LENGTH) + ' and '
                          + str(REF_SEQ_MAX_LENGTH) + ' nucleotides.'))

    return Sequence(ref_sequence)

def validate_custom_primers(cfp: str, crp: str, ref_sequence, hsps):
    """ Verify that the custom primers contain a valid alphabet. If
    both are provided, then return then as Primer objects. Else,
    verify that the one that is provided has exactly one binding site
    in the reference sequence and none in the non-target regions.
    """
    if not isinstance(cfp, (str, type(None))):
        raise TypeError('Custom forward primer must be a string or None.')

    if not isinstance(crp, (str, type(None))):
        raise TypeError('Custom reverse primer must be a string or None.')

    if cfp is None:
        cfp = ''

    if crp is None:
        crp = ''

    cfp = Sequence(cfp.upper())
    crp = Sequence(crp.upper())

    if regex.search('[^ACGT]', str(cfp)):
        raise ValueError(('Custom forward primer contains invalid characters. '
                          'It must only contain the characters A, C, G, or '
                          'T.'))

    if regex.search('[^ACGT]', str(crp)):
        raise ValueError(('Custom forward primer contains invalid characters. '
                          'It must only contain the characters A, C, G, or '
                          'T.'))

    if cfp and crp:
        # If both are provided, then we just calculated the PCR conditions
        # of these two. No reference sequence is required.
        cfp = Primer(cfp, 0, 0, 1)
        crp = Primer(crp, 0, 0, -1)
        return cfp, crp

    hsp_seqs = list(map(lambda hsp: hsp.s_seq, hsps))

    if cfp:
        sites = (binding_sites((ref_sequence,), cfp)
                 + binding_sites((ref_sequence,), cfp.rev_comp()))
        if len(sites) == 0:
            raise NestedLoopError(('Forward primer has no binding site in the '
                                   'reference sequence.'))
        elif len(sites) > 1:
            raise NestedLoopError(('Forward primer has multiple binding sites '
                                   'in the reference sequence.'))

        hsp_sites = (binding_sites(hsp_seqs, cfp)
                     + binding_sites(hsp_seqs, cfp.rev_comp()))
        if len(hsp_sites) > 0:
            raise NestedLoopError(('Forward primer has binding sites in the '
                                   'non-target regions.'))
        cfp = Primer(sequence=cfp, span=(sites[0].start(), sites[0].end()),
                     strand=1, custom=True)
    else:
        cfp = ''

    if crp:
        sites = (binding_sites((ref_sequence,), crp)
                 + binding_sites((ref_sequence,), crp.rev_comp()))
        if len(sites) == 0:
            raise NestedLoopError(('Reverse primer has no binding site in the '
                                   'reference sequence.'))
        elif len(sites) > 1:
            raise NestedLoopError(('Reverse primer has multiple binding sites '
                                   'in the reference sequence.'))

        hsp_sites = (binding_sites(hsp_seqs, crp)
                     + binding_sites(hsp_seqs, crp.rev_comp()))
        if len(hsp_sites) > 0:
            raise NestedLoopError(('Reverse primer has binding sites in the '
                                   'non-target regions.'))
        crp = Primer(sequence=crp, span=(sites[0].start(), sites[0].end()),
                     strand=-1, custom=True)
    else:
        crp = ''

    return cfp, crp

def validate_amplicon_size(pcr_min: int, pcr_max: int):
    from nestedloop import PCR_MINIMUM, PCR_MAXIMUM

    if not isinstance(pcr_min, int) or not isinstance(pcr_max, int):
        raise TypeError('Amplicon sizes must be defined by integers.')

    if (not PCR_MINIMUM <= pcr_min <= PCR_MAXIMUM
            or not PCR_MINIMUM <= pcr_max <= PCR_MAXIMUM):
        raise ValueError(('Desired amplicon size must be between '
                          + str(PCR_MINIMUM) + ' and ' + str(PCR_MAXIMUM)
                          + ' base pairs.'))

    if pcr_min >= pcr_max:
        raise ValueError(('Minimum amplicon size must be smaller than the '
                          'maximum amplicon size.'))

    return pcr_min, pcr_max

def validate_fspan(f_from: int, f_to: int, ref_sequence):
    from nestedloop import F_FROM_MIN, F_TO_MIN

    if not isinstance(f_from, int) or not isinstance(f_to, int):
        raise TypeError('Regions must be defined by integers.')

    if not F_FROM_MIN <= f_from or not F_TO_MIN <= f_to:
        raise ValueError('Forward primer regions must be positive.')

    if f_from >= f_to:
        raise ValueError('\'From\' indices must be smaller than \'to\' indices.')

    f_to = len(ref_sequence) if f_to > len(ref_sequence) else f_to

    # Zero-index the region.
    return f_from-1, f_to

def validate_rspan(r_from: int, r_to: int, ref_sequence):
    from nestedloop import R_FROM_MIN, R_TO_MIN

    if not isinstance(r_from, int) or not isinstance(r_to, int):
        raise TypeError('Regions should be defined by integers.')

    if not R_FROM_MIN <= r_from or not R_TO_MIN <= r_to:
        raise ValueError('Reverse primer regions must be positive.')

    if r_from >= r_to:
        raise ValueError('\'From\' indices must be smaller than \'to\' indices.')

    r_to = len(ref_sequence) if r_to > len(ref_sequence) else r_to

    # Zero-index the region.
    return r_from-1, r_to

def validate_temperature(tm: float):
    from nestedloop import TEMPERATURE_MIN, TEMPERATURE_MAX
    if not isinstance(tm, (int, float)):
        raise TypeError('Temperatures should be numeric.')

    if not TEMPERATURE_MIN <= tm <= TEMPERATURE_MAX:
        raise ValueError(('Temperature values should be between '
                          + str(TEMPERATURE_MIN) + ' and '
                          + str(TEMPERATURE_MAX) + ' degrees celsius.'))

    return tm

def validate_num_to_return(num_to_return: int):
    if not isinstance(num_to_return, int):
        raise TypeError('Only an integer number of pairs can be returned.')

    if num_to_return < 1:
        raise ValueError('The number of pairs to return must be positive.')

    return num_to_return

def validate_nontargets(nontargets: str):
    if not isinstance(nontargets, (str, type(None))):
        raise TypeError('Non-targets must be a string or None.')

    if nontargets:
        if nontargets.startswith('<?xml version="1.0"?>'):
            parser = XmlParser()
        else:
            parser = PairwiseParser()
        hsps = parser.parse(nontargets)
    else:
        hsps = list()

    return hsps
