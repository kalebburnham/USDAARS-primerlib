def validate_input_data(data: str) -> str:
    """ Validates and returns a cleaned string of the user's input data
    with SNPs. Note that the SNP grammar is validated later.
    """
    from starp import DATA_MAX_LENGTH

    if not isinstance(data, str):
        raise TypeError('Input must be a string.')

    # Remove FASTA header.
    if data.startswith('>'):
        data = ''.join(data.splitlines()[1:])

    data = data.strip().upper()

    if len(data) > DATA_MAX_LENGTH:
        raise ValueError(f'Input data must be less than {DATA_MAX_LENGTH} '
                          'nucleotides.')

    return data
