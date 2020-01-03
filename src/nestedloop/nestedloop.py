"""
License information goes here.
"""

import re

from .models import Primer, Sequence
from .exceptions import NestedLoopError
from .parsers import XmlParser, PairwiseParser
from .utils import binding_sites, combine, primer_generate, primer_sort

TEMPERATURE_MIN = 35.0  # In Degrees Celsius
TEMPERATURE_MAX = 72.5
NUM_TO_RETURN_MIN = 1
PCR_MINIMUM = 1  # Minimum Amplicon Size
PCR_MAXIMUM = 20000  # Maximum Amplicon Size
REF_SEQ_MIN_LENGTH = 100
REF_SEQ_MAX_LENGTH = 20000
F_FROM_MIN = 1
F_TO_MIN = 1
R_FROM_MIN = 1
R_TO_MIN = 1

class NestedLoop:
    """ The only public class for NestedLoop. This is the main file,
    and everything gets called from here.

    Attributes:
        ref_sequence - A string of the reference sequence, possibly in FASTA format.
        tm - A 3-tuple of (minimum, optimum, maximum) temperature.
        f_from - An int representing the first point to look for f primers.
        f_to - An int representing the last point to look for f primers.
        r_from - An int representing the first point to look for r primers.
        r_to - An int representing the last point to look for r primers.
        pcr_min - An integer, the minimum distance between the primers.
        pcr_max - An integer, the maximum distance between the primers.
        num_to_return - An integer, the maximum number of pairs to return, if possible.
        custom_forward_primer - A string of any custom f primer. Else, empty string.
            Should be written 5'->3' on PLUS strand.
        custom_reverse_primer - A string of any custom r primer. Else, empty string.
            Should be written 5'->3' on MINUS strand.

    """

    def __init__(self, ref_sequence, tm, f_from, f_to, r_from, r_to,
                 pcr_min, pcr_max, num_to_return, custom_forward_primer=None,
                 custom_reverse_primer=None, nontargets=None):
        """
        Arguments:
            See above for the class attributes.
            nontargets: The BLAST alignment data.
            ** Regions should be entered into the program as 1-indexed inclusive
        """
        self.ref_sequence = validate_reference_sequence(ref_sequence)
        self.tm_min = validate_temperature(tm[0])
        self.tm_opt = validate_temperature(tm[1])
        self.tm_max = validate_temperature(tm[2])
        self.f_from, self.f_to = validate_fspan(f_from, f_to, ref_sequence)
        self.r_from, self.r_to = validate_rspan(r_from, r_to, ref_sequence)
        self.pcr_min, self.pcr_max = validate_amplicon_size(pcr_min, pcr_max)
        self.num_to_return = validate_num_to_return(num_to_return)
        self.hsps = validate_nontargets(nontargets) # High-Scoring Pairs
        self.custom_forward_primer, self.custom_reverse_primer = (
                        validate_custom_primers(custom_forward_primer,
                                                custom_reverse_primer,
                                                self.ref_sequence,
                                                self.hsps))

        self.f_primers = list()
        self.r_primers = list()
        self.pairs = list()

    def run(self):

        if self.custom_forward_primer:
            self.f_primers = [self.custom_forward_primer]
        else:
            self.f_primers = primer_generate(self.ref_sequence, 1, self.f_from, self.f_to, (self.tm_min, self.tm_opt, self.tm_max))

        if self.custom_reverse_primer:
            self.r_primers = [self.custom_reverse_primer]
        else:
            self.r_primers = primer_generate(ref_sequence=self.ref_sequence, strand=-1,
                                             start=self.r_from, stop=self.r_to,
                                             tms=(self.tm_min, self.tm_opt, self.tm_max))

        self.f_primers = primer_sort(self.f_primers)
        self.r_primers = primer_sort(self.r_primers)

        self.pairs = combine(self.f_primers, self.r_primers,
                             self.num_to_return, self.pcr_min, self.pcr_max,
                             self.ref_sequence, self.hsps)

        for pair in self.pairs:
            pair.additive(self.ref_sequence)

        return self.pairs

def validate_reference_sequence(ref_sequence: str):
    """ Removes the FASTA header (if there is one) from the input which is
    assumed to be on the first line then strips any new line characters or
    whitespace at the beginning and end. """
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

    if re.search('[^ACGT]', str(cfp)):
        raise ValueError(('Custom forward primer contains invalid characters. '
                          'It must only contain the characters A, C, G, or '
                          'T.'))

    if re.search('[^ACGT]', str(crp)):
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
        cfp = Primer(cfp, sites[0].start(), sites[0].end(), strand=1, custom=True)
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
        crp = Primer(crp, sites[0].start(), sites[0].end(), strand=-1, custom=True)
    else:
        crp = ''

    return cfp, crp

def validate_amplicon_size(pcr_min: int, pcr_max: int):
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
