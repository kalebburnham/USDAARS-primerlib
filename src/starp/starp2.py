"""
License information goes here.
"""

import textwrap
import html

from .amasfactory import AmasFactory
from .utils import (rgenerate, rfilter, rfilter_complementary, rsorted,
    segregate, rfilter_by_binding_sites, rtailed)
from .models import Sequence, Snp
from .parsers import get_parser
from .exceptions import StarpError

class Starp2:
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
        self.input_data = input_data.upper()

        self.hsps = []
        self.amas = None
        self.reverse_primers = []
        self.num_to_return = 5
        self.pcr_max = 500 ### ********* TEMPORARY VALUE

        self.snp_parser = get_parser(self.input_data)
        self.snp = None  # The user's chosen Snp.
        self.snps = self.snp_parser.snps()

        # The lengths of each allele may be different. They only
        # contain nucleotides or 'N'.
        self.allele1 = self.snp_parser.allele1
        self.allele2 = self.snp_parser.allele2

        # The lengths of each aligned sequence are equal. These may
        # contain dashes to represent insertions and deletions.
        self.allele1_aligned = Sequence(self.snp_parser.allele1_aligned)
        self.allele2_aligned = Sequence(self.snp_parser.allele2_aligned)

    def choose_snp(self, snp):
        """
        Given a user's chosen Snp, automatically proceeds to generate
        AMAS and reverse primers.

        Args:
            snp: A Snp object representing the user's chosen SNP.

        Returns:
            None.

        Raises:
            Value Error: the chosen SNP is not recognized.
            Type Error: the passed argument snp is not a Snp object.
            Starp Error: the snp is a substitution, and there is
                another SNP within 4-25 bp upstream.
        """
        if not snp in self.snps:
            raise ValueError("Chosen SNP is not in the known list of SNPs.")

        if not isinstance(snp, Snp):
            raise TypeError("Chosen SNP is not of correct type.")

        if snp.type == 'substitution':
            nearby_snps = [s for s in self.snps
                           if snp.start-25 <= s.start <= snp.start-25]
            for nearby_snp in nearby_snps:
                if nearby_snp.type != "substitution":
                    raise StarpError("""Cannot create AMAS primers here due to
                                        the presence of insertions/deletions
                                        in the region.""")

        self.snp = snp

        # Create the two amas primers. They are stored in a 2-tuple.
        amas_factory = AmasFactory(self.allele1, self.allele2, self.snps, snp)
        self.amas = amas_factory.generate()

        # Create reverse primers.
        candidates = rgenerate(self.allele1, snp,
                                  min_length=18, max_length=27)

        candidates = rfilter(candidates, self.amas, self.snps, self.pcr_max)

        low_tm_primers, high_tm_primers = segregate(candidates)

        # Binding sites will be checked on these sequences.
        sequences = (self.allele1, self.allele2) + tuple(self.hsps)

        if high_tm_primers:
            high_tm_primers = rfilter_complementary(high_tm_primers)
            high_tm_primers = rsorted(high_tm_primers)

            
            self.reverse_primers = rfilter_by_binding_sites(high_tm_primers,
                                                       self.allele1,
                                                       self.allele2,
                                                       self.hsps,
                                                       max_num=3,
                                                       amas=self.amas)

        # This branch is reached when the high_tm_primers does not contain
        # any acceptable primers.
        if low_tm_primers and not self.reverse_primers:
            tailed_rprimers = rtailed(low_tm_primers)
            breakpoint()
            tailed_rprimers = rfilter_complementary(tailed_rprimers)
            breakpoint()
            tailed_rprimers = rsorted(tailed_rprimers)
            breakpoint()
            self.reverse_primers = rfilter_by_binding_sites(low_tm_primers,
                                                            self.allele1,
                                                            self.allele2,
                                                            self.hsps,
                                                            max_num=3,
                                                            amas=self.amas)
            breakpoint()

        if not self.reverse_primers:
            raise StarpError("No Reverse Primers Found")

        return self.reverse_primers


    def html(self):
        """
        Converts the two alleles to a human-readable, blast-like HTML
        format with SNPs converted to HTML buttons.

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

        allele1_lines = textwrap.wrap(str(self.allele1_aligned), WIDTH, break_on_hyphens=False)
        allele2_lines = textwrap.wrap(str(self.allele2_aligned), WIDTH, break_on_hyphens=False)
        midline_lines = textwrap.wrap(midline, WIDTH, replace_whitespace=False)

        lines_in_order = []

        for i in range(len(allele1_lines)):
            lines_in_order.append(allele1_lines[i])
            lines_in_order.append(midline_lines[i].replace('-', ' '))
            lines_in_order.append(allele2_lines[i])
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
