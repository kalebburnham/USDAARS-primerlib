"""
License information goes here.
"""

import textwrap
import html
import logging

from .amasfactory import (generate_amas_for_substitution,
                          generate_amas_for_indel, substitute_bases)
from .utils import (rgenerate, rfilter, rfilter_complementary, rsorted,
                    segregate, rfilter_by_binding_sites, rtailed, add_tails,
                    record_spans)
from .models import Sequence, Snp, AmasGroup
from .parsers import get_parser
from .exceptions import StarpError
from .data_validation import validate_input_data

__all__ = [
    'amasfactory', 'exceptions', 'models', 'parsers', 'utils',
]

DATA_MAX_LENGTH = 10000

class Starp:
    """
    Attributes:
        input_data: The user's input.
        hsps:
        amas:
        reverse_primers:
        num_to_return: The number of reverse primers to return
        snp: The user's chosen SNP, around which AMAS primers will be
            created.
        snps: A list of all SNPs between the two alleles.
        allele1: The allele with all the first choices of each SNP
        allele2: The allele with all the second choices of each SNP
    """
    def __init__(self, input_data):
        """
        Initializes Starp data.

        Args:
            input_data: The raw user data with SNPs.

        Returns:
            None

        Raises:
            None

        """
        self.input_data = validate_input_data(input_data)

        self.hsps = []
        self.amas = None
        self.upstream = None
        self.downstream = None
        self.reverse_primers = []
        self.num_to_return = 5
        self.pcr_max = 500 ### ********* TEMPORARY VALUE

        self.snp = None  # The user's chosen Snp object.

        snp_parser = get_parser(self.input_data)
        self.snps = snp_parser.snps()

        # The lengths of each allele may be different. They only
        # contain nucleotides or 'N'.
        self.allele1 = snp_parser.allele1
        self.allele2 = snp_parser.allele2

        # The lengths of each aligned sequence are equal. These may
        # contain dashes to represent insertions and deletions.
        self.allele1_aligned = Sequence(snp_parser.allele1_aligned)
        self.allele2_aligned = Sequence(snp_parser.allele2_aligned)

        logging.getLogger().setLevel(logging.INFO)

    def set_snp(self, descriptor: str):
        self.snp = Snp(descriptor)

    def run(self):
        if not self.snp:
            raise ValueError(('Chosen snp has not been set. Use '
                              'starp.set_snp(descriptor) and try again.'))

        if not isinstance(self.snp, Snp):
            raise TypeError("Chosen SNP is not of correct type.")

        if self.snp not in self.snps:
            raise ValueError("Chosen SNP is not in the known list of SNPs.")

        logging.info('Generating AMAS primers.')

        if self.snp.type == 'deletion':
            # In deletion snps, the descriptor does not contain information on
            # what the reference nucleotide is, so it must be added..
            self.snp.ref_nucleotide = str(self.allele1[self.snp.position])

        if self.snp.type == 'substitution':
            upstream_amas, downstream_amas = generate_amas_for_substitution(
                self.allele1, self.allele2,
                self.snp.position, self.snps
            )
            upstream_amas = substitute_bases(upstream_amas, self.snp,
                                             'upstream')
            downstream_amas = substitute_bases(downstream_amas, self.snp,
                                               'downstream')

        elif self.snp.type == 'insertion' or self.snp.type == 'deletion':
            upstream_amas, downstream_amas = generate_amas_for_indel(
                self.allele1, self.allele2,
                self.snp.position, self.snps
            )
            upstream_amas = substitute_bases(upstream_amas, self.snp,
                                             'upstream')
            downstream_amas = substitute_bases(downstream_amas, self.snp,
                                               'downstream')

        logging.info('Finished generating AMAS primers.')
        logging.info('Making reverse primers.')

        # Create reverse primers.
        candidates = rgenerate(self.allele1,
                               min_length=18, max_length=27)

        # Downstream candidates are those AFTER the snp, so they pair
        # with the upstream amas pair, and vice versa.

        upstream_rcandidates = rfilter(candidates, downstream_amas,
                                       self.snps, self.pcr_max)
        downstream_rcandidates = rfilter(candidates, upstream_amas,
                                         self.snps, self.pcr_max)

        upstream = AmasGroup(upstream_amas, downstream_rcandidates)
        upstream.rprimers = record_spans(upstream.rprimers, self.allele1,
                                         self.allele2)
        upstream.add_rtails()
        upstream.rprimers = rsorted(upstream.rprimers)
        upstream.rprimers = rfilter_by_binding_sites(upstream.rprimers,
                                                     self.allele1,
                                                     self.allele2,
                                                     self.hsps,
                                                     max_num=3,
                                                     amas=upstream.amas)

        downstream = AmasGroup(downstream_amas, upstream_rcandidates)
        downstream.rprimers = record_spans(downstream.rprimers, self.allele1,
                                           self.allele2)
        downstream.add_rtails()
        downstream.rprimers = rsorted(downstream.rprimers)
        downstream.rprimers = rfilter_by_binding_sites(downstream.rprimers,
                                                       self.allele1,
                                                       self.allele2,
                                                       self.hsps,
                                                       max_num=3,
                                                       amas=downstream.amas)

        if not upstream.rprimers and not downstream.rprimers:
            raise StarpError('No reverse primers found.')

        self.downstream = downstream
        self.upstream = upstream


        # START HERE
        #self.amas = add_tails()

    def html(self):
        """
        Converts the two alleles to a human-readable, blast-like HTML
        format with SNPs converted to HTML buttons.

        Output will look like this:

        GTACTGATGACTGACTGCGCGCGCGCCGGGGGGGCTGAGATGCAGTCGTACGTCAAGCAA
        ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        GTACTGATGACTGACTGCGCGCGCGCCGGGGGGGCTGAGATGCAGTCGTACGTCAAGCAA

        GTACTGATGACTGACTGCGCGCGCGCCGGGGGGGCTGAGATGCAGTCGTACGTCAAGCAA
        |||||||||||||||||||||||||||||||||||||||||||||  |||||||||||||
        GTACTGATGACTGACTGCGCGCGCGCCGGGGGGGCTGAGATGCAGAAGTACGTCAAGCAA

        Args:
            None

        Returns:
            A preformatted HTML string that resembles BLAST output.
            Each SNP is converted to a button.

        Raises:
            None
        """
        WIDTH = 60

        midline = ''.join(['|' if a == b else '-'
                           for a, b in zip(str(self.allele1_aligned),
                                           str(self.allele2_aligned))])

        allele1_lines = textwrap.wrap(str(self.allele1_aligned),
                                      width=WIDTH,
                                      break_on_hyphens=False)
        allele2_lines = textwrap.wrap(str(self.allele2_aligned),
                                      width=WIDTH,
                                      break_on_hyphens=False)
        midline_lines = textwrap.wrap(midline,
                                      width=WIDTH,
                                      replace_whitespace=False)

        lines_in_order = []

        for line_num, line in enumerate(allele1_lines):
            lines_in_order.append(allele1_lines[line_num])
            lines_in_order.append(midline_lines[line_num].replace('-', ' '))
            lines_in_order.append(allele2_lines[line_num])
            lines_in_order.append('')

        button_style = '"font-size: 1em; padding: 0; \
                         border: none; cursor: pointer; background-color: lightgreen;"'
        alignment = '<pre><p>' + '\n'.join(lines_in_order) + '</p></pre>'
        # Add buttons to the string. Assumes self.snps are in order.
        snps = iter(self.snps)
        alignment = list(alignment)

        for idx, char in enumerate(alignment, 0):
            if char == ' ':
                snp = next(snps)
                alignment[idx] = ('<button class="snp" name="snp" value='
                                  + html.escape(snp.descriptor)
                                  + ' style='+button_style+'> </button>')

        alignment = ''.join(alignment)

        return alignment
