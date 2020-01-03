"""
License information goes here.
"""

from .exceptions import StarpError
from .models import Sequence, Primer

class AmasFactory:
    """
    A factory for creating AMAS primers.

    Attributes:
        allele1: The continuous allele 1 Sequence from the SNP
                parsers.
        allele2: The continuous allele 2 Sequence from the SNP
                parsers.
        snps: The list of all SNPs between these alleles.
        snp: The user's chosen SNP.

    Public methods:
        generate

    """

    def __init__(self, allele1, allele2, snps, snp):
        """

        Args:
            allele1: The continuous allele 1 Sequence from the SNP
                parsers. Should not contain any dashes.
            allele2: The continuous allele 2 Sequence from the SNP
                parsers. Should not contain any dashes.
            snps: The list of all SNPs between these alleles.
            snp: The user's chosen SNP.

        Raises:
            None
        """
        self.allele1 = Sequence(allele1)
        self.allele2 = Sequence(allele2)
        self.snps = snps
        self.snp = snp

    def generate(self):
        """
        Generates an AMAS pair using rules defined by Dr. Long.

        Args:
            None

        Returns:
            A 2-tuple holding the AMAS Sequence objects.

        Raises:
            StarpError: when no suitable Starp primers are found.
            StarpError: when a non-substitution SNP is selected.
        """
        amas = (None, None)

        if self.snp.type == 'substitution':
            # This information came from
            # docs/Long_11_Nov_19.docx
            pairs = self._generate_upstream(16, 26)
            pair = _best_pair(pairs)

            if pair is None:
                pairs = self._generate_downstream(16, 26)
                pair = _best_pair(pairs)

            if pair is None:
                raise StarpError('Cannot find Starp primers at this location '
                                 'due to improper melting temperatures.')

            amas = self._substitute_bases(pair)

        elif self.snp.type == 'insertion' or self.snp.type == 'deletion':
            # The information from this branch came from
            # docs/how to design AMAS primers for indel_20191125[3981].pptx
            pairs = self._generate_downstream(0, 8)

            # Remove pairs with same base at 3' end.
            pairs = filter(lambda pair: pair[0][-1] != pair[1][-1], pairs)

            # If pairs is empty after the filter, try upstream.
            if not pairs:
                pairs = self._generate_upstream(0, 8)
                pairs = filter(lambda pair: pair[0][-1] != pair[1][-1], pairs)

            # Order the pairs based on
            # 1) The number of nucleotide differences in the last 4 bases.
            #    More differences is preferable.
            # 2) Length of the sequences in each pair. Longer is preferable.
            pairs = sorted(pairs,
                           key=lambda pair: (
                               Sequence.hamming(pair[0][-4:], pair[1][-4:]),
                               len(pair[0])),
                           reverse=True)

            amas = pairs[0]

        else:
            raise StarpError('AMAS creation is undefined for non-substitution '
                             'SNPs.')

        return amas

    def _generate_upstream(self, minimum, maximum):
        """
        Generates F primer pairs upstream from the SNP site.
        For example, if minimum = 16 and maximum = 25, we grab the 16
        through 25 nucleotides BEFORE the SNP.
        These would be the pairs that get created.

        F1                NNNNNNNNNNNNNNNNA
        F2                NNNNNNNNNNNNNNNNB

        F1               NNNNNNNNNNNNNNNNNA
        F2               NNNNNNNNNNNNNNNNNB

        F1              NNNNNNNNNNNNNNNNNNA
        F2              NNNNNNNNNNNNNNNNNNB

        F1             NNNNNNNNNNNNNNNNNNNA
        F2             NNNNNNNNNNNNNNNNNNNB

        F1            NNNNNNNNNNNNNNNNNNNNA
        F2            NNNNNNNNNNNNNNNNNNNNB

        F1           NNNNNNNNNNNNNNNNNNNNNA
        F2           NNNNNNNNNNNNNNNNNNNNNB

        F1          NNNNNNNNNNNNNNNNNNNNNNA
        F2          NNNNNNNNNNNNNNNNNNNNNNB

        F1         NNNNNNNNNNNNNNNNNNNNNNNA
        F2         NNNNNNNNNNNNNNNNNNNNNNNB

        F1        NNNNNNNNNNNNNNNNNNNNNNNNA
        F2        NNNNNNNNNNNNNNNNNNNNNNNNB

        F1       NNNNNNNNNNNNNNNNNNNNNNNNNA
        F2       NNNNNNNNNNNNNNNNNNNNNNNNNB

        In the case of insertions or deletions, see
            'docs/how to design AMAS-primers for Indel*'

        Args:
            minimum: The shortest distance upstream to traverse.
            maximum: The longest distance upstream to traverse.

        Returns:
            A list of tuples (F1, F2).

        Raises:
            None.
        """

        # Find the first index on the upstream that is not a '-'
        # in allele1.
        pos1 = self.snp.start-1
        while self.allele1[pos1] == '-':
            pos1 -= 1

        # Find the first index on the upstream that is not a '-'
        # in allele2.
        pos2 = self.snp.start-1
        while self.allele2[pos2] == '-':
            pos2 -= 1

        # Negative indices are allowed in Python, but are no good for
        # primer creation.
        if pos1 < 0 or pos2 < 0:
            raise StarpError("Cannot create AMAS primers at this position.")

        # Create the two lists and zip them together.
        allele1_primers = [Primer(self.allele1[pos1-i:pos1+1], pos1-i, pos1+1, 1)
                             for i in range(minimum, maximum)]
        allele2_primers = [Primer(self.allele2[pos2-i:pos2+1], pos2-i, pos2+1, 1) 
                             for i in range(minimum, maximum)]
        pairs = zip(allele1_primers, allele2_primers)
        print(bool(pairs))
        return list(pairs)

    def _generate_downstream(self, minimum, maximum):
        """
        Generates F primer pairs downstream from the SNP site. The
        shortest primer length is minimum+1 and the longest is
        maximum+1.

        If minimum=16 and maximum=25, These are the pairs that get
        created:

        F1       ANNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNNNNNN

        F1       ANNNNNNNNNNNNNNNNNNNNNNNNN
        F2       BNNNNNNNNNNNNNNNNNNNNNNNNN


        Args:
            None. Uses the self attributes of the class.

        Returns:
            A list of tuples (F1, F2).

        """

        # Find the first position on the downstream that is not a '-'
        # character in allele1.
        pos1 = self.snp.start-1
        while self.allele1[pos1] == '-':
            pos1 += 1

        # Find the first position on the downstream that is not a '-'
        # character in allele2.
        pos2 = self.snp.start-1
        while self.allele2[pos2] == '-':
            pos2 += 1

        # Negative indices are allowed in Python, but are no good for
        # primer creation.
        if pos1 < 0 or pos2 < 0:
            raise StarpError("Cannot create AMAS primers at this position.")

        # Create the two lists and zip them together.
        allele1_sequences = [Primer(self.allele1[pos1:pos1+i+1], pos1, pos1+i+1, 1)
                             for i in range(minimum, maximum)]
        allele2_sequences = [Primer(self.allele2[pos2:pos2+i+1], pos2, pos2+i+1, 1)
                             for i in range(minimum, maximum)]
        pairs = zip(allele1_sequences, allele2_sequences)

        return list(pairs)

    def _substitute_bases(self, pair):
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
        the substitution index may be gleaned from the matrix.

        Args:
            pair: A tuple with the selected AMAS primer pair.

        Returns:
            A tuple holding the AMAS pair Sequences with substituted
            bases.

        Raises:
            StarpError: when a non-substitution SNP is selected.
            StarpError: when a non-substitution SNP is within last 4
                nucleotides of the 3' end.
        """

        # Other snps in the 2nd, 3rd, or 4th positions from the chosen snp.
        extra_snps = [s for s in self.snps
                      if s.start > self.snp.start-4
                      and s.start < self.snp.start]

        snp = self.snp
        amas1 = pair[0]
        amas1_str = str(pair[0])
        amas2 = pair[1]
        amas2_str = str(pair[1])
        if snp.type != "substitution":
            raise StarpError("Cannot substitute bases on a non-substitution snp.")

        if not extra_snps:
            # There are no extra SNPs closeby to the chosen one.

            sub_index = {frozenset(('C', 'G')) :
                             {'GGGC' : -3, 'GGGG' : -2, 'GGAC' : (-3, -4), 'GGAG' : (-3, -4),
                              'GAGC' : (-2, -4), 'GAGG' : (-2, -4), 'GAAC' : (-2, -3), 'GAAG' : (-2, -3),
                              'AGGC' : (-2, -3), 'AGGG' : (-2, -3), 'AGAC' : (-2, -4), 'AGAG' : (-2, -4),
                              'AAGC' : (-3, -4), 'AAGG' : (-3, -4), 'AAAC' : -3, 'AAAG' : -2},
                         frozenset(('C', 'T')) :
                             {'GGGC' : -2, 'GGGT' : -3, 'GGAC' : -4, 'GGAT' : -3,
                              'GAGC' : -2, 'GAGT' : -4, 'GAAC' : -2, 'GAAT' : -3,
                              'AGGC' : -2, 'AGGT' : -3, 'AGAC' : -2, 'AGAT' : -4,
                              'AAGC' : -4, 'AAGT' : -3, 'AAAC' : -2, 'AAAT' : -3},
                         frozenset(('C', 'A')) :
                             {'GGGC' : -2, 'GGGA' : -3, 'GGAC' : -3, 'GGAA' : -4,
                              'GAGC' : -2, 'GAGA' : -4, 'GAAC' : -2, 'GAAA' : -3,
                              'AGGC' : -2, 'AGGA' : -3, 'AGAC' : -2, 'AGAA' : -4,
                              'AAGC' : -3, 'AAGA' : -4, 'AAAC' : -2, 'AAAA' : -3},
                         frozenset(('G', 'T')) :
                             {'GGGG' : -2, 'GGGT' : -3, 'GGAG' : -4, 'GGAT' : -3,
                              'GAGG' : -2, 'GAGT' : -4, 'GAAG' : -2, 'GAAT' : -4,
                              'AGGG' : -2, 'AGGT' : -3, 'AGAG' : -2, 'AGAT' : -4,
                              'AAGG' : -4, 'AAGT' : -3, 'AAAG' : -2, 'AAAT' : -2},
                         frozenset(('G', 'A')) :
                             {'GGGG' : -2, 'GGGA' : -3, 'GGAG' : -3, 'GGAA' : -4,
                              'GAGG' : -2, 'GAGA' : -4, 'GAAG' : -2, 'GAAA' : -3,
                              'AGGG' : -2, 'AGGA' : -3, 'AGAG' : -2, 'AGAA' : -4,
                              'AAGG' : -3, 'AAGA' : -4, 'AAAG' : -2, 'AAAA' : -3},
                         frozenset(('A', 'T')) :
                             {'GGGT' : -3, 'GGGA' : -4, 'GGAT' : -3, 'GGAA' : -4,
                              'GAGT' : -2, 'GAGA' : -4, 'GAAT' : -2, 'GAAA' : -3,
                              'AGGT' : -2, 'AGGA' : -3, 'AGAT' : -2, 'AGAA' : -4,
                              'AAGT' : -3, 'AAGA' : -4, 'AAAT' : -3, 'AAAA' : -4}}

            # There are some patterns to Dr. Long's substitution instructions.
            #
            # 1) He only considers the last 4 nucleotides.
            # 2) He regularly refers to G and C together, and A at T together. So,
            #    mapping these to a common value reduces the number of keys we
            #    need to consider. The mapping is: C,G -> G and A,T -> T.
            # 3) The allele also matters. In a G/C SNP, one allele has a G and
            #    one has a C, for example. When generate_upstream() is called,
            #    this value is naturally at the end. However, generate_downstream()
            #    has the value at the start. For consistency, place it at the end.
            #    This reduces the size of the dictionary by half.

            # This key structure encodes a lot of information. His two pages of
            # solid text for the case of no other SNPs at the 3' end can be
            # defined in just 29 (long) lines of code.

            allele1_key = (amas1_str[-4:-1].replace('C', 'G').replace('T', 'A')
                           + snp.ref_nucleotide)
            allele2_key = (amas2_str[-4:-1].replace('C', 'G').replace('T', 'A')
                           + snp.new_nucleotide)

            amas1 = Primer(_substitute(amas1_str, sub_index[snp.nucleotides][allele1_key]),
                           amas1.start, amas1.end, amas1.strand)
            amas2 = Primer(_substitute(amas2_str, sub_index[snp.nucleotides][allele2_key]),
                           amas2.start, amas2.end, amas2.strand)

        elif len(extra_snps) == 1:
            other = extra_snps[0]
            if other.type is not "substitution":
                raise StarpError("Cannot substitute bases if a non-substitution snp is nearby.")

            # What the 4-tuple means:
            # Let it be ('CG', -4, 'AT', -3)
            # Then, if the allele has a C or G in the additional SNP position,
            # then substitute the -4th index. If it has an A or T in the
            # additional SNP position, substitute the -3 index.

            sub_index = {frozenset(('C', 'G')) :
                             {'GGPC' : -4, 'GGPG' : -3, 'GGNC' : -4, 'GGNG' : -3,
                              'GAPC' : -4, 'GAPG' : -3, 'GANC' : ('CG', -4, 'AT', -3), 'GANG' : ('CG', -4, 'AT', -3),
                              'AGPC' : -4, 'AGPG' : -3, 'AGNC' : ('CG', -3, 'AT', -4), 'AGNG' : ('CG', -3, 'AT', -4),
                              'AAPC' : -4, 'AAPG' : -3, 'AANC' : -4, 'AANG' : -3,
                              'GPGC' : -4, 'GPGG' : -2, 'GPAC' : -4, 'GPAG' : -2,
                              'GNGC' : -4, 'GNGG' : -2, 'GNAC' : ('CG', -4, 'AT', -2), 'GNAG' : ('CG', -4, 'AT', -2),
                              'APGC' : -4, 'APGG' : -2, 'APAC' : -4, 'APAG' : -2,
                              'ANGC' : ('CG', -2, 'AT', -4), 'ANGG' : ('CG', -2, 'AT', -4), 'ANAC' : -4, 'ANAG' : -2,
                              'PGGC' : -3, 'PGGG' : -2, 'PGAC' : -3, 'PGAG' : -2,
                              'PAGC' : -3, 'PAGG' : -2, 'PAAC' : -3, 'PAAG' : -2,
                              'NGGC' : -3, 'NGGG' : -2, 'NGAC' : ('CG', -3, 'AT', -2), 'NGAG' : ('CG', -3, 'AT', -2),
                              'NAGC' : ('CG', -2, 'AT', -3), 'NAGG' : ('CG', -2, 'AT', -3), 'NAAC' : -3, 'NAAG' : -2},
                         frozenset(('C', 'T')) : {
                             'GGPC' : -4, 'GGPT' : -3, 'GGNC' : -4, 'GGNT' : -3,
                             'GAPC' : -4, 'GAPT' : -3, 'GANC' : ('CG', -4, 'AT', -4), 'GANT' : ('CG', -3, 'AT', -3),
                             'AGPC' : -3, 'AGPT' : -4, 'AGNC' : ('CG', -3, 'AT', -4), 'AGNT' : ('CG', -3, 'AT', -4),
                             'AAPC' : -4, 'AAPT' : -3, 'AANC' : -4, 'AANT' : -3,
                             'GPGC' : -2, 'GPGT' : -4, 'GPAC' : -4, 'GPAT' : -2,
                             'GNGC' : -2, 'GNGT' : -4, 'GNAC' : ('CG', -4, 'AT', -2), 'GNAT' : ('CG', -4, 'AT', -2),
                             'APGC' : -2, 'APGT' : -4, 'APAC' : -2, 'APAT' : -4,
                             'ANGC' : ('CG', -2, 'AT', -2), 'ANGT' : ('CG', -4, 'AT', -4), 'ANAC' : -2, 'ANAT' : -4,
                             'PGGC' : -2, 'PGGT' : -4, 'PGAC' : -3, 'PGAT' : -2,
                             'PAGC' : -2, 'PAGT' : -3, 'PAAC' : -2, 'PAAT' : -3,
                             'NGGC' : -2, 'NGGT' : -3, 'NGAC' : ('CG', -3, 'AT', -2), 'NGAT' : ('CG', -3, 'AT', -2),
                             'NAGC' : ('CG', -2, 'AT', -2), 'NAGT' : ('CG', -3, 'AT', -3), 'NAAC' : -2, 'NAAT' : -3},
                         frozenset(('C', 'A')) : {
                             'GGPC' : -3, 'GGPA' : -4, 'GGNC' : -3, 'GGNA' : -4,
                             'GAPC' : -4, 'GAPA' : -3, 'GANC' : ('CG', -4, 'AT', -3), 'GANA' : ('CG', -4, 'AT', -3),
                             'AGPC' : -3, 'AGPA' : -4, 'AGNC' : ('CG', -3, 'AT', -3), 'AGNA' : ('CG', -4, 'AT', -4),
                             'AAPC' : -3, 'AAPA' : -4, 'AANC' : -3, 'AANA' : -4,
                             'GPGC' : -2, 'GPGA' : -4, 'GPAC' : -4, 'GPAA' : -2,
                             'GNGC' : -2, 'GNGA' : -4, 'GNAC' : ('CG', -4, 'AT', -2), 'GNAA' : ('CG', -4, 'AT', -2),
                             'APGC' : -2, 'APGA' : -4, 'APAC' : -2, 'APAA' : -4,
                             'ANGC' : ('CG', -2, 'AT', -2), 'ANGA' : ('CG', -4, 'AT', -4), 'ANAC' : -2, 'ANAA' : -4,
                             'PGGC' : -2, 'PGGA' : -3, 'PGAC' : -3, 'PGAA' : -2,
                             'PAGC' : -2, 'PAGA' : -3, 'PAAC' : -2, 'PAAA' : -3,
                             'NGGC' : -2, 'NGGA' : -3, 'NGAC' : ('CG', -3, 'AT', -2), 'NGAA' : ('CG', -3, 'AT', -2),
                             'NAGC' : ('CG', -2, 'AT', -2), 'NAGA' : ('CG', -3, 'AT', -3), 'NAAC' : -2, 'NAAA' : -3},
                         frozenset(('G', 'T')) : {
                             'GGPT' : -3, 'GGPG' : -4, 'GGNT' : -3, 'GGNG' : -4,
                             'GAPT' : -3, 'GAPG' : -4, 'GANT' : ('CG', -3, 'AT', -3), 'GANG' : ('CG', -4, 'AT', -4),
                             'AGPT' : -4, 'AGPG' : -3, 'AGNT' : ('CG', -3, 'AT', -4), 'AGNG' : ('CG', -3, 'AT', -4),
                             'AAPT' : -3, 'AAPG' : -4, 'AANT' : -3, 'AANG' : -4,
                             'GPGT' : -4, 'GPGG' : -2, 'GPAT' : -2, 'GPAG' : -4,
                             'GNGT' : -4, 'GNGG' : -2, 'GNAT' : ('CG', -4, 'AT', -2), 'GNAG' : ('CG', -4, 'AT', -2),
                             'APGT' : -4, 'APGG' : -2, 'APAT' : -4, 'APAG' : -2,
                             'ANGT' : ('CG', -4, 'AT', -4), 'ANGG' : ('CG', -2, 'AT', -2), 'ANAT' : -4, 'ANAG' : -2,
                             'PGGT' : -3, 'PGGG' : -2, 'PGAT' : -2, 'PGAG' : -3,
                             'PAGT' : -3, 'PAGG' : -2, 'PAAT' : -3, 'PAAG' : -2,
                             'NGGT' : -3, 'NGGG' : -2, 'NGAT' : ('CG', -3, 'AT', -2), 'NGAG' : ('CG', -3, 'AT', -2),
                             'NAGT' : ('CG', -3, 'AT', -3), 'NAGG' : ('CG', -2, 'AT', -2), 'NAAT' : -3, 'NAAG' : -2},
                         frozenset(('G', 'A')) : {
                             'GGPA' : -4, 'GGPG' : -3, 'GGNA' : -4, 'GGNG' : -3,
                             'GAPA' : -3, 'GAPG' : -4, 'GANA' : ('CG', -4, 'AT', -3), 'GANG' : ('CG', -4, 'AT', -3),
                             'AGPA' : -4, 'AGPG' : -3, 'AGNA' : ('CG', -4, 'AT', -4), 'AGNG' : ('CG', -3, 'AT', -3),
                             'AAPA' : -4, 'AAPG' : -3, 'AANA' : -4, 'AANG' : -3,
                             'GPGA' : -4, 'GPGG' : -2, 'GPAA' : -2, 'GPAG' : -4,
                             'GNGA' : -4, 'GNGG' : -2, 'GNAA' : ('CG', -4, 'AT', -2), 'GNAG' : ('CG', -4, 'AT', -2),
                             'APGA' : -4, 'APGG' : -2, 'APAA' : -4, 'APAG' : -2,
                             'ANGA' : ('CG', -4, 'AT', -4), 'ANGG' : ('CG', -2, 'AT', -2), 'ANAA' : -4, 'ANAG' : -2,
                             'PGGA' : -3, 'PGGG' : -2, 'PGAA' : -2, 'PGAG' : -3,
                             'PAGA' : -3, 'PAGG' : -2, 'PAAA' : -3, 'PAAG' : -2,
                             'NGGA' : -3, 'NGGG' : -2, 'NGAA' : ('CG', -3, 'AT', -2), 'NGAG' : ('CG', -3, 'AT', -2),
                             'NAGA' : ('CG', -3, 'AT', -3), 'NAGG' : ('CG', -2, 'AT', -2), 'NAAA' : -3, 'NAAG' : -2},
                         frozenset(('A', 'T')) : {
                             'GGPA' : -4, 'GGPT' : -3, 'GGNA' : -4, 'GGNT' : -3,
                             'GAPA' : -4, 'GAPT' : -3, 'GANA' : ('CG', -4, 'AT', -3), 'GANT' : ('CG', -4, 'AT', -3),
                             'AGPA' : -4, 'AGPT' : -3, 'AGNA' : ('CG', -3, 'AT', -4), 'AGNT' : ('CG', -3, 'AT', -4),
                             'AAPA' : -4, 'AAPT' : -3, 'AANA' : -4, 'AANT' : -3,
                             'GPGA' : -4, 'GPGT' : -2, 'GPAA' : -4, 'GPAT' : -2,
                             'GNGA' : -4, 'GNGT' : -2, 'GNAA' : ('CG', -4, 'AT', -2), 'GNAT' : ('CG', -4, 'AT', -2),
                             'APGA' : -4, 'APGT' : -2, 'APAA' : -4, 'APAT' : -2,
                             'ANGA' : ('CG', -2, 'AT', -4), 'ANGT' : ('CG', -2, 'AT', -4), 'ANAA' : -4, 'ANAT' : -2,
                             'PGGA' : -3, 'PGGT' : -2, 'PGAA' : -3, 'PGAT' : -2,
                             'PAGA' : -3, 'PAGT' : -2, 'PAAA' : -3, 'PAAT' : -2,
                             'NGGA' : -3, 'NGGT' : -2, 'NGAA' : ('CG', -3, 'AT', -2), 'NGAT' : ('CG', -3, 'AT', -2),
                             'NAGA' : ('CG', -2, 'AT', -3), 'NAGT' : ('CG', -2, 'AT', -2), 'NAAA' : -3, 'NAAT' : -2}}

            # See the case above for an explanation of this.
            allele1_key = (amas1_str[-4:-1].replace('C', 'G').replace('T', 'A')
                           + snp.ref_nucleotide)
            allele2_key = (amas2_str[-4:-1].replace('C', 'G').replace('T', 'A')
                           + snp.new_nucleotide)

            # Other SNP relative position.
            # If the other SNP is at the 3rd index from the end, then this
            # value is set to -3, for example.
            other_snp_rel_pos = other.start - snp.start - 1

            # Convert the allele*_end to a format which can query the
            # substitution dictionary.
            if other.nucleotides == {'C', 'G'} or other.nucleotides == {'A', 'T'}:
                placeholder = 'P'
            else:
                placeholder = 'N'

            allele1_key = list(allele1_key)
            allele1_key[other_snp_rel_pos] = placeholder
            allele1_key = ''.join(allele1_key)

            allele2_key = list(allele2_key)
            allele2_key[other_snp_rel_pos] = placeholder
            allele2_key = ''.join(allele2_key)

            index_to_substitute1 = sub_index[allele1_key]
            index_to_substitute2 = sub_index[allele2_key]

            # Deal with the output of sub_index if it returns a tuple
            # like ('CG', -2, 'AT', -3).
            if isinstance(index_to_substitute1, tuple):
                sub_iter = iter(index_to_substitute1)
                if other.ref_nucleotide == next(sub_iter):
                    index_to_substitute1 = next(sub_iter)
                else:
                    next(sub_iter)
                    next(sub_iter)
                    index_to_substitute1 = next(sub_iter)

            if isinstance(index_to_substitute2, tuple):
                sub_iter = iter(index_to_substitute2)
                if other.new_nucleotide == next(sub_iter):
                    index_to_substitute2 = next(sub_iter)
                else:
                    next(sub_iter)
                    next(sub_iter)
                    index_to_substitute2 = next(sub_iter)

            amas1 = Primer(_substitute(amas1_str, index_to_substitute1), amas1.start, amas1.end, amas1.strand)
            amas2 = Primer(_substitute(amas2_str, index_to_substitute2), amas2.start, amas2.end, amas2.strand)

        else:
            # No substitutions necessary.
            pass

        return (amas1, amas2)

def _best_pair(pairs):
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

def _substitute(allele, idx):
    """
    Substitutes the given index of the allele according to the
    mapping nucleotide_sub. This is necessary because strings
    and Sequence objects are immutable.

    Args:
        allele: the Sequence object to modify
        idx: an index or tuple of the indices to change.

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
