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
    """ Return the 'hseqs' (hit sequences) from the nontarget xml data.
    """
    nontarget_seqs = []

    if nontargets:
        if nontargets.startswith('<?xml'):
            parser = XmlParser()
        else:
            raise ValueError('Nontarget data not in XML format.')
            # parser = PairwiseParser()
        nontarget_seqs = parser.parse(nontargets)

    return nontarget_seqs
