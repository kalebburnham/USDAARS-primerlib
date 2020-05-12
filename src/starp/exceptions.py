class StarpError(Exception):
    """ A basic error to distinguish between program issues and
    conditions that forces Starp to stop. """
    def __init__(self, message):
        self.message = message
