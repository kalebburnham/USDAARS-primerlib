from .sequence import Sequence
from .primer import Primer

class AmasPrimer(object):
    """
    A container for an AMAS primer.

    Public Methods:
        melting_temperature()


     Represents a normal primer with a tail. 
    Question: Should length, Tm, and GC calculations be done on the
    original primer? In that case, should this inherit from Primer? """

    def __init__(self, primer, tail):
        """ Tail should just be an integer, 1 or 2. """
        self.primer = primer
        self.tail = tail

    def __len__(self):
        return len(self.amas)

    def __str__(self):
        return self.amas.sequence

    def melting_temperature(self):
        """
        Gives the melting temperature of the original
        """
        return self.primer.melting_temperature()

    def gc(self):
        return self.primer.gc

    def html(self):
        """ One document shows the tail underlined and the substitutions
        highlighted. """

