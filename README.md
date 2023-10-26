# StrucSeq
Fetch and process data between protein structures and sequences.

## Functions

### get_uniprot_accessions

get_uniprot_accessions(pdbcode, strict, selenium (non functional), debug)

Takes a pdb assession code and accesses the PDB to get the Uniprot assession code of the protein for each chain

#### Parameters

pdbcode : str  
The 4-letter protein data bank (PDB) code.

strict : bool, optional  
If False, will not raise an exception in input is incorrect.

selenium : bool, optional  
Choose to use selenium to scrabe the browser version, rather than the XML version (currently not got this feature working)

debug : bool, optional  
Write information to console

#### Returns

dict  
Uniprot accession codes of each chain in the strcture. Key is the chain identifier,
content is the accession code.

#### Examples 

    get_uniprot_accessions("3ii6")
    
    {'A': 'Q13426',
    'B': 'Q13426',
    'C': 'Q13426',
    'D': 'Q13426',
    'X': 'P49917',
    'Y': 'P49917'}

### iterate_uniprot_accessions

iterate_uniprot_accessions(in_csv, chain_cols, out_csv, delimiter, debug):

Takes an input CSV with a PDBid header and specified custom header(s) to get the
uniprot ID of the chains and saves a CSV with the assession codes and unique structures
and chains.

#### Parameters

in_csv : str  
The the input CSV with PDBids that is used to generate accession codes.

chain_cols : list or str  
Column or list of columns that should be used to identify the chain(s) when
fetching their accessions.  
e.g. ["a chain", "b chain"]

out_csv : str  
The output CSV used to store fetched accession codes.

delimiter : str, optional  
Delimiter of the CSV file, default "\t" tab.

debug : bool, optional  
Whether to print progess or other notifications

#### Returns

Nothing is returned, but a CSV is saved specifying the uniprot ID of each chain in every structure of the input CSV.