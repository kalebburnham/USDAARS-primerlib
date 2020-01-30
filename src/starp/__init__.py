"""
License information goes here.
"""

import textwrap
import html
import logging
import itertools
import time

from .amasfactory import (generate_amas_for_substitution,
                          generate_amas_for_indel, substitute_bases)
from .utils import (rgenerate, rfilter, rsorted, rfilter_by_binding_sites,
                    rtailed, add_tails, aligned_position)
from .models import Sequence, Snp, AmasGroup
from .parsers import get_parser
from .exceptions import StarpError
from .data_validation import validate_input_data, validate_nontargets

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
    def __init__(self, input_data, nontargets=None):
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
        
        hsps = validate_nontargets(nontargets)
        self.nontargets = [hsp.s_seq for hsp in hsps]

        self.amas = None
        self.upstream_pairs = list()
        self.downstream_pairs = list()
        self.reverse_primers = list()
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

        snp_index = aligned_position(self.snp.position, self.snps)

        logging.info('Generating AMAS primers.')

        if self.snp.type == 'deletion':
            # In deletion snps, the descriptor does not contain information on
            # what the reference nucleotide is, so it must be added..
            self.snp.ref_nucleotide = str(self.allele1[self.snp.position])

        # Generate the best AMAS pair depending on the type of polymorphism.
        if self.snp.type == 'substitution':
            upstream_amas, downstream_amas = generate_amas_for_substitution(
                self.allele1_aligned, self.allele2_aligned,
                self.snp.position
            )

        elif self.snp.type == 'insertion' or self.snp.type == 'deletion':
            upstream_amas, downstream_amas = generate_amas_for_indel(
                self.allele1_aligned, self.allele2_aligned,
                self.snp.position
            )

        # Substitute bases according to
        # docs/starp/STARP F primer design[4312].docx pages 3-13.
        upstream_amas = substitute_bases(upstream_amas, self.snp,
                                         snp_position='last')
        downstream_amas = substitute_bases(downstream_amas, self.snp,
                                           snp_position='first')

        if not upstream_amas and not downstream_amas:
            raise StarpError('Could not generate AMAS primers.')

        logging.info('Finished generating AMAS primers.')

        logging.info('Making reverse primers.')

        # Create reverse primers common to both alleles. This does not
        # guarantee unique binding sites.
        candidates = rgenerate(self.allele1_aligned, self.allele2_aligned,
                               min_length=18, max_length=27)

        # Downstream candidates are those AFTER the snp, so they pair
        # with the upstream amas pair, and vice versa.

        # *********
        # Attempt to pair upstream AMAS primers with downstream rprimers.

        # I couldn't think of very good names for the containers holding
        # the amas primers and the rprimers, so here is an explanation.

        # amas_r_group is the container holding AMAS primers upstream from
        # the SNP and rprimers downstream of the SNP, while r_amas_group is
        # the container holding rprimers upstream from the SNP and amas
        # primers downstream from the SNP. 

        # This gives a little intuition if you imagine in amas_r_group that
        # the amas primers come before the r primers, and vice versa for
        # r_amas_group.
        if upstream_amas:
            downstream_rcandidates = rfilter(candidates, upstream_amas,
                                             self.pcr_max, snp_position='last')

            amas_r_group = AmasGroup(upstream_amas, downstream_rcandidates)

            start_time = time.time()
            amas_r_group.rprimers = rfilter_by_binding_sites(amas_r_group.rprimers,
                                                         self.allele1,
                                                         self.allele2,
                                                         self.nontargets,
                                                         amas=amas_r_group.amas)
            logging.info((f'Filtered upstream rprimers by binding sites in '
                          f'{str(time.time()-start_time)} seconds.'))
            amas_r_group.add_rtails()
            amas_r_group.rprimers = rsorted(amas_r_group.rprimers)[:3]
        else:
            # Need to define something, else errors arise because variable
            # does not exist.
            amas_r_group = AmasGroup([], [])

        # List of amas primers with common reverse primer.
        # [(amas1, amas2), rprimer]
        if amas_r_group:
            amas_r_group.rprimers = [primer if primer.strand == -1 else primer.rev_comp()
                                     for primer in amas_r_group.rprimers]

            # Cartesian product of the two sets. The AMAS primers need to be
            # paired with each individual rprimer to add on the correct tail.
            amas_r_pairs = list(itertools.product([amas_r_group.amas],
                                                    amas_r_group.rprimers))
            for pair in amas_r_pairs:
                amplicon1 = abs(pair[0][0].start - pair[1].allele1_end)
                amplicon2 = abs(pair[0][1].start - pair[1].allele2_end)
                add_tails(*pair[0], amplicon1, amplicon2, snp_position='last')

            self.upstream_pairs = amas_r_pairs
        else:
            self.upstream_pairs = []

        # ********
        # Attempt to pair downstream amas with upstream rprimers.

        if downstream_amas:
            upstream_rcandidates = rfilter(candidates, downstream_amas,
                                           self.pcr_max, snp_position='first')

            r_amas_group = AmasGroup(downstream_amas, upstream_rcandidates)

            start_time = time.time()
            r_amas_group.rprimers = rfilter_by_binding_sites(r_amas_group.rprimers,
                                                           self.allele1,
                                                           self.allele2,
                                                           self.nontargets,
                                                           amas=r_amas_group.amas)
            logging.info((f'Filtered downstream rprimers by binding sites in '
                          f'{str(time.time()-start_time)} seconds.'))
            r_amas_group.add_rtails()
            r_amas_group.rprimers = rsorted(r_amas_group.rprimers)[:3]
            
        else:
            r_amas_group = AmasGroup([], [])

        if r_amas_group:
            r_amas_group.rprimers = [primer if primer.strand == -1 else primer.rev_comp()
                                   for primer in r_amas_group.rprimers]

            r_amas_pairs = list(itertools.product([r_amas_group.amas],
                                                      r_amas_group.rprimers))
            for pair in r_amas_pairs:
                amplicon1 = abs(pair[0][0].end - pair[1].allele1_start)
                amplicon2 = abs(pair[0][1].end - pair[1].allele2_start)
                add_tails(*pair[0], amplicon1, amplicon2, snp_position='first')

            self.downstream_pairs = r_amas_pairs
        else:
            self.downstream_pairs = []

        if not amas_r_group.rprimers and not r_amas_group.rprimers:
            raise StarpError('No reverse primers found.')

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
