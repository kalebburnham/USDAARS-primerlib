import re

from .models import Primer, Sequence
from .exceptions import NestedLoopError
from .parsers import XmlParser, PairwiseParser
from .utils import binding_sites, combine, primer_generate, primer_sort
from .data_validation import (validate_reference_sequence,
                              validate_custom_primers,
                              validate_amplicon_size,
                              validate_fspan,
                              validate_rspan,
                              validate_temperature,
                              validate_num_to_return,
                              validate_nontargets)

__all__ = [
    'additive', 'exceptions', 'nestedloop', 'parsers', 'models'
]

TEMPERATURE_MIN = 35.0  # Degrees Celsius
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
