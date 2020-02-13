from .parsers import XmlParser, PairwiseParser

def validate_input_data(data: str) -> str:
    """ Validates and returns a cleaned string of the user's input data
    with SNPs. Note that the SNP grammar is validated later.
    """
    from starp import DATA_MAX_LENGTH

    if not isinstance(data, str):
        raise TypeError('Input must be a string.')

    data = data.strip().upper()

    if len(data) > DATA_MAX_LENGTH:
        raise ValueError(f'Input data must be less than {DATA_MAX_LENGTH} '
                          'characters.')

    return data

def validate_nontargets(nontargets: str):
    """ Attempts to parse user input alignment information and
    returns HSP objects. """
    if not isinstance(nontargets, (str, type(None))):
        raise TypeError('Non-targets must be a string or None.')

    if nontargets:
        if nontargets.startswith('<?xml version="1.0"?>'):
            parser = XmlParser()
        else:
            parser = PairwiseParser()
        hsps = parser.parse(nontargets)
    else:
        hsps = list()

    return hsps
