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
        
        self.nontargets = validate_nontargets(nontargets)

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
        snp = Snp(descriptor)
        if snp.type == 'deletion':
            snp.ref_nucleotide = str(self.allele1_aligned[snp.position])
        self.snp = snp

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
        upstream_amas = substitute_bases(upstream_amas, snp_position='last')
        downstream_amas = substitute_bases(downstream_amas, snp_position='first')

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
            # Since rprimers come after the SNP, they need to be
            # reverse complemented to place them on the minus strand.
            amas_r_group.rprimers = [primer if primer.strand == -1 else primer.rev_comp()
                                     for primer in amas_r_group.rprimers]

            # Cartesian product of the two sets. The AMAS primers need to be
            # paired with each individual rprimer to add on the correct tail.
            amas_r_pairs = list(itertools.product([amas_r_group.amas],
                                                    amas_r_group.rprimers))
            for pair in amas_r_pairs:
                # Amplicon length on first allele.
                amplicon1 = abs(pair[0][0].start - pair[1].allele1_end)

                # Amplicon length on second allele.
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
            # Format of r_amas_group:
            # [[[amas1, amas2], [rprimer1, rprimer2, primer3, ...]]

            # Format of r_amas_pairs:
            # [[[amas1, amas2], rprimer], [[amas1, amas2], rprimer2], ...]
            r_amas_pairs = list()

            for pair in list(itertools.product([r_amas_group.amas], r_amas_group.rprimers)):
                # Amplicon length on first allele.
                amplicon1 = abs(pair[0][0].end - pair[1].allele1_start)

                # Amplicon length on second allele.
                amplicon2 = abs(pair[0][1].end - pair[1].allele2_start)
                add_tails(*pair[0], amplicon1, amplicon2, snp_position='first')

                # Since AMAS primers come after the SNP, they need to be
                # reverse complemented to place them on the minus strand.
                r_amas_pairs.append([[pair[0][0].rev_comp(), pair[0][1].rev_comp()], pair[1]])

            self.downstream_pairs = r_amas_pairs
        else:
            self.downstream_pairs = []

        if not amas_r_group.rprimers and not r_amas_group.rprimers:
            raise StarpError('No reverse primers found.')

    def html(self, width=60):
        """
        Converts the two alleles to a human-readable, blast-like HTML
        format with SNPs converted to HTML buttons.

        Output will look like this:

        GTACTGATGACTGACTGCNCGCGCGCCGGGGGGGCTGAGATGCAGTCGTACGTCAAGCAA
        |||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||
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

        # Create the midline where '|' signifies that nucleotides are
        # the same, 'N' means one of the nucleotides is an 'N', and '-'
        # signifies a SNP. The dashes get replaced to spaces later, but
        # textwrap.wrap does not work well with whitespace.
        midline = []
        for a, b in zip(str(self.allele1_aligned), str(self.allele2_aligned)):
            if a == 'N' or b == 'N':
                char = '*'
            elif a == b:
                char = '|'
            else:
                char = '-'

            midline.append(char)

        midline = ''.join(midline)

        # Wrap the lines to the specified width.
        allele1_lines = textwrap.wrap(str(self.allele1_aligned),
                                      width=width,
                                      break_on_hyphens=False)
        allele2_lines = textwrap.wrap(str(self.allele2_aligned),
                                      width=width,
                                      break_on_hyphens=False)
        midline_lines = textwrap.wrap(midline,
                                      width=width,
                                      replace_whitespace=False)

        # Create an array of the lines that resembles the output. Eg,
        # line 0:   GTACTGTGC
        # line 1:   |||| ||*|
        # line 2:   GTACAGTNC
        # line 3:   <empty line>
        lines_in_order = []

        for line_num, line in enumerate(allele1_lines):
            lines_in_order.append(allele1_lines[line_num])
            lines_in_order.append(midline_lines[line_num].replace('-', ' '))
            lines_in_order.append(allele2_lines[line_num])
            lines_in_order.append('')

        # The HTML style of the buttons.
        button_style = '"font-size: 1em; padding: 0; \
                         border: none; cursor: pointer; background-color: lightgreen;"'
        
        # Signify that this is preformatted HTML code so it is not
        # escaped in the browser. 'alignment' is the entire HTML block
        # that will be returned.
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

        # Replace any 'N's with whitespace. This must be done at the
        # end so it is not mistaken for a SNP.
        alignment = alignment.replace('*', ' ')

        return alignment
