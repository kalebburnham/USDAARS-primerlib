"""
License information goes here.
"""

import math

class Additive:
    """
    A class that automatically calculates PCR conditions given a
    target sequence upon initialization.

    Refer to docs/NL-PCR_reaction_conditions.docx

    Attributes:
        sequence: The user-defined target sequence.
        max_gc: The maximum gc percentage of any 100bp segment
        min_gc: The minimum gc percentage of any 100bp segment.
        additive: The amount of additive needed in PCR.
        pcr_temperatures: The temperature sequence for PCR.
        n: Unknown. Refer to docs/NL-PCR_reaction_conditions.docx
    """

    def __init__(self, target_sequence):
        self.sequence = target_sequence
        self.max_gc = self._max_gc()
        self.min_gc = self._min_gc()
        self.additive = self._additive()
        self.pcr_temperatures = self._pcr_temperatures()

        # What is n? It is on the PCR conditions document.
        self.n = 1 + math.ceil((len(target_sequence) / 3000) - 0.3)

    def _additive(self):
        """ Units are in ul. """
        additive = None

        if math.isnan(self.max_gc):
            return float('nan')

        if self.max_gc >= 0.95:
            additive = "6.5-9.0"
        elif self.max_gc >= 0.90:
            additive = "6.0-8.5"
        elif self.max_gc >= 0.85:
            additive = "5.5-8.0"
        elif self.max_gc >= 0.80:
            additive = "3.5-6.0"
        elif self.max_gc >= 0.75:
            additive = "3.0-5.5"
        elif self.max_gc >= 0:
            additive = "2.5-5.0"

        return additive

    def _max_gc(self):
        """ Returns the max GC content of any 100 bp subsequence
        in sequence. """

        # What to do about Ns?
        if len(self.sequence) < 100:
            # What to do?
            return float('nan')

        # Get all 100 bp sequence from the target sequence.
        fragments = [self.sequence[i:i+100] for i in range(len(self.sequence)-99)]
        gcs = [fragment.gc for fragment in fragments]
        return max(gcs)

    def _min_gc(self):
        if len(self.sequence) < 100:
            return float('nan')

        fragments = [self.sequence[i:i+100] for i in range(len(self.sequence)-99)]
        gcs = [fragment.gc for fragment in fragments]
        return min(gcs)

    def _pcr_temperatures(self):
        """ Returns temperatures in a pre-formatted HTML form."""
        temperatures = None

        if self.max_gc >= 0.95:
            if self.min_gc <= 0.10:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 58/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 59/30\"\n"
            else:
                temperatures = "1) 72/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
        elif self.max_gc >= 0.90:
            if self.min_gc <= 0.10:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 59/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
            else:
                temperatures = "1) 71/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 61/30\"\n"
        elif self.max_gc >= 0.85:
            if self.min_gc <= 0.10:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 60/30\"\n"
            elif self.min_gc <= 0.30:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n6) 61/30\"\n"
            else:
                temperatures = "1) 70/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n"
        elif self.max_gc >= 0.80:
            temperatures = "1) 69/30\"\n2) 68/30\"\n3) 66/30\"\n4) 64/30\"\n5) 62/30\"\n"
        else:
            temperatures = "1) 68/30\"\n2) 66/30\"\n3) 64/30\"\n4) 62/30\"\n"
            #temperatures = "68/30\"; 66/30\"; 64/30\"; 62/30\";"

        return temperatures
