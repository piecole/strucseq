# StrucSeq
Fetch and process data between protein structures and sequences.

## Functions

get_uniprot_accessions(pdbcode, strict, selenium (non functional), debug)

Takes a pdb assession code and accesses the PDB to get the Uniprot assession code of the protein for each chain

    Parameters

    pdbcode : str
        The 4-letter protein data bank (PDB) code.
    strict : bool, optional
        If False, will not raise an exception in input is incorrect.
    selenium : bool, optional
        Choose to use selenium to scrabe the browser version, rather than the XML 
        version (currently not got this feature working)
    debug : bool, optional
        Write information to console

    Returns

    dict
        Uniprot accession codes of each chain in the strcture. Key is the chain identifier,
        content is the accession code.

    Example
    <pre><code>get_uniprot_accessions("3ii6")<br>
    {'A': 'Q13426',
    'B': 'Q13426',
    'C': 'Q13426',
    'D': 'Q13426',
    'X': 'P49917',
    'Y': 'P49917'}
    </code></pre>