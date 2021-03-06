import textwrap
import html
import itertools
import time

from .amasfactory import (generate_amas_for_substitution,
                          generate_amas_for_indel)
from .utils import (add_rtails, rgenerate, rfilter, rfilter_tailed_primers,
                    rsorted, rfilter_by_binding_sites,
                    rtailed, add_tails)
from .models import Sequence, Snp, StarpGroup
from .parsers import get_parser
from .data_validation import validate_input_data, validate_nontargets

__all__ = [
    'amasfactory', 'exceptions', 'models', 'parsers', 'utils',
]

DATA_MAX_LENGTH = 10000

class Starp:
    """
    Attributes:
        input_data: The user's input.
        nontargets: A list of hit sequences. Could be empty.
        snp: The user's chosen SNP, around which AMAS primers will be
            created.
        snps: A list of all SNPs between the two alleles.
        allele1: The allele with all the first choices of each SNP
        allele2: The allele with all the second choices of each SNP
        allele1_aligned: The alignment of allele1 with insertions and
                         deletions.
        allele2_aligned: The alignment of allele2 with insertions and
                         deletions.
        starp_groups: 
    """
    def __init__(self, input_data, nontargets=None):
        """
        Initializes Starp data.

        Args:
            input_data: The raw user data with SNPs.

        Returns:
            None
        """
        self.input_data = validate_input_data(input_data)
        
        self.nontargets = validate_nontargets(nontargets)

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

        self.starp_groups = []

        # Generate a Starp Group for each SNP.
        # A StarpGroup contains a SNP, an AMAS primer for both alleles,
        # and the relative SNP position to both alleles (first or last
        # nucleotide). Later, the groups receive rcandidates and each
        # group determines which are appropriate rprimers for itself.
        # These are saved in self.starp_groups
        self.prototype()

    def set_snp(self, descriptor: str):
        snp = Snp(descriptor)
        if snp.type == 'deletion':
            snp.ref_nucleotide = str(self.allele1_aligned[snp.position])
        self.snp = snp

    def prototype(self):
        """ Generate Starp Groups for each SNP and add them to
        self.starp_groups. """
        for snp in self.snps:
            # Generate the best AMAS pair depending on the type of polymorphism.
            if snp.type == 'substitution':
                upstream_amas, downstream_amas = generate_amas_for_substitution(
                    self.allele1_aligned, self.allele2_aligned,
                    snp.position
                )
            elif snp.type == 'insertion' or snp.type == 'deletion':
                upstream_amas, downstream_amas = generate_amas_for_indel(
                    self.allele1_aligned, self.allele2_aligned,
                    snp.position
                )

            if upstream_amas:
                self.starp_groups.append(StarpGroup(upstream_amas[0], upstream_amas[1],
                                         snp_position='last', snp=snp))

            if downstream_amas:
                self.starp_groups.append(StarpGroup(downstream_amas[0], downstream_amas[1],
                                         snp_position='first', snp=snp))

        return self.starp_groups

    def run(self):
        # Create reverse primers common to both alleles. This does not
        # guarantee unique binding sites.
        rcandidates = rgenerate(self.allele1_aligned, self.allele2_aligned,
                               min_length=18, max_length=27)

        # Filter the primers based on alphabet, repetitive nucleotides,
        # and very high/low GC content.
        rcandidates = rfilter(rcandidates)

        # Remove rcandidates with multiple binding sites.
        rcandidates = rfilter_by_binding_sites(rcandidates, self.allele1,
                                               self.allele2, self.nontargets)

        # Add tails to low melting temperature rprimers.
        rcandidates = add_rtails(rcandidates)

        # Remove rprimers that became too long, have too high of a melting
        # temperature, or significant complementary scores after modifications
        # from add_rtails.
        rcandidates = rfilter_tailed_primers(rcandidates)

        # Sort rprimers from best to worst.
        rcandidates = rsorted(rcandidates)

        tail1 = Sequence('GCAACAGGAACCAGCTATGAC')
        tail2 = Sequence('GACGCAAGTGAGCAGTATGAC')

        # A list of StarpGroups with tails added to AMAS primers.
        tailed_groups = []

        # The following steps are completed for each starp group:
        # 1. The group is supplied with all possible rcandidates.
        # 2. Sort the rcandidates and select the best primers for
        #    this group's AMAS primer pair.
        # 3. Split the group into 2 groups. In group1, the AMAS primers
        #    are attached tail1 and tail2 respectively. In group2, the
        #    AMAS primers are attached tail2 and tail1.
        #    
        #    Then, each group keeps the rprimers that agree with the
        #    tails attached to the group's AMAS primers.
        #
        # 4. If each group has suitable primers, add it to tailed_groups.
        for group in self.starp_groups:
            group.rcandidates = rcandidates
            group.set_rprimers()
            group1, group2 = group.segregate(tail1, tail2)

            if group1.rprimers:
                tailed_groups.append(group1)
            if group2.rprimers:
                tailed_groups.append(group2)

        # Primers that come before the SNP should be on the plus strand,
        # while those after the primer should be on the minus strand.
        for group in tailed_groups:
            if group.snp_position == 'first':
                # The AMAS primers need to be rev comped to be placed
                # on the -1 strand.
                if group.amas1.strand == 1:
                    group.amas1 = group.amas1.rev_comp()
                if group.amas2.strand == 1:
                    group.amas2 = group.amas2.rev_comp()
            elif group.snp_position == 'last':
                # The reverse primers need to be rev comped.
                group.rprimers = [primer.rev_comp() for primer in group.rprimers]

        # Replace Starp groups with the groups with AMAS tails.
        self.starp_groups = tailed_groups
        return self.starp_groups

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
