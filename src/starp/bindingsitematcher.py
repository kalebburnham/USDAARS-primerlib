"""
License information goes here.
"""

import regex

class BindingSiteMatcher:
    """
    Given a reference sequence and another sequence, attempts to find
    all binding sites of the second sequence in the reference sequence.

    A 'binding site' is defined as a place where there are <= tolerance
    mismatches. Both the ref_sequence and its reverse complement are
    searched.
    """

    def match(self, ref_sequence, seq, tolerance=0, stop=10):
        """
        Args:
            ref_sequence: The reference Sequence.
            seq: A smaller Sequence to search for in ref_sequence.
            tolerance: The number of acceptable nucleotide mismatches for a match.
            stop: Stop searching if this many matches are found.
                Defaults to 10.

        Returns:
            A list of starp.Match objects.

        Raises:
            None
        """
        binding_sites = []

        pattern = '(' + str(seq) + '){s<=' + str(tolerance) + '}'
        rev_comp_pattern = '(' + str(seq.rev_comp()) + '){s<=' + str(tolerance) + '}'

        match_iter = regex.finditer(pattern,
                                    str(ref_sequence),
                                    overlapped=True)

        rev_comp_match_iter = regex.finditer(rev_comp_pattern,
                                             str(ref_sequence),
                                             overlapped=True)

        for match in match_iter:
            if len(binding_sites) >= stop:
                break
            binding_sites.append(Match(seq, match.group(), match.span()))


        for match in rev_comp_match_iter:
            if len(binding_sites) >= stop:
                break
            binding_sites.append(Match(seq.rev_comp(), match.group(), match.span()))

        return binding_sites

class Match:
    """
    A container for matches found by BindingSiteMatcher.

    Attributes:
        seq: The original seq that was searched for.
        match: The string sequence that was matched against the
            ref_sequence.
        span: A tuple containing the beginning and ending indices in
            the binding site match.

    Public Methods:
        hamming
        positions
    """

    def __init__(self, seq, match, span):
        self.seq = seq
        self.match = match
        self.span = span

    def __eq__(self, other):
        if not isinstance(other, Match):
            return False

        return (self.seq == other.seq
                and self.match == other.match
                and self.span == other.span)

    def hamming(self):
        """
        Returns the hamming distance between seq and match. This is the
        number of mismatches between them.
        """
        assert len(self.seq) == len(self.match)
        return sum(n1 != n2 for n1, n2 in zip(self.seq, self.match))

    def positions(self):
        """
        Returns a list of booleans with True whenever seq and match
        contain the same indices.
        """
        return [n1 == n2 for n1, n2 in zip(self.seq, self.match)]
        