class NestedLoopError(Exception):
    """
    Generic error that will immediately stop the regular execution
    of the Nested Loop program.
    """
    def __init__(self, message):
        self.message = message

class BindingSiteError(NestedLoopError):
    def __init__(self, message):
        self.message = message
