# StrucSeq
A python package to fetch and process data between protein structures and sequences.

This is very much a work in progress.

# Functions

## get_uniprot_accessions

```python
get_uniprot_accessions(pdbcode, strict, selenium (non functional), debug)
```

Takes a pdb assession code and accesses the PDB to get the Uniprot assession code of the protein for each chain

### Parameters

**pdbcode** : str  
The 4-letter protein data bank (PDB) code.

**strict** : bool, optional  
If False, will not raise an exception in input is incorrect.

**selenium** : bool, optional  
Choose to use selenium to scrabe the browser version, rather than the XML version (currently not got this feature working)

**debug** : bool, optional  
Write information to console

### Returns

**dict**  
Uniprot accession codes of each chain in the strcture. Key is the chain identifier,
content is the accession code.

### Examples

```python
>>> get_uniprot_accessions("3ii6")

{'A': 'Q13426',
'B': 'Q13426',
'C': 'Q13426',
'D': 'Q13426',
'X': 'P49917',
'Y': 'P49917'}
```

## iterate_uniprot_accessions

```python
iterate_uniprot_accessions(in_csv, chain_cols, out_csv, delimiter, debug):
```

Takes an input CSV with a PDBid header and specified custom header(s) to get the
uniprot ID of the chains and saves a CSV with the assession codes and unique structures
and chains.

### Parameters

**in_csv** : str  
The the input CSV with PDBids that is used to generate accession codes.

**chain_cols** : list or str  
Column or list of columns that should be used to identify the chain(s) when
fetching their accessions.  
e.g. ["a chain", "b chain"]

**out_csv** : str  
The output CSV used to store fetched accession codes.

**delimiter** : str, optional  
Delimiter of the CSV file, default "\t" tab.

**debug** : bool, optional  
Whether to print progess or other notifications

### Returns

Nothing is returned, but a CSV is saved specifying the uniprot ID of each chain in every structure of the input CSV.


## get_uniprot_details

```python
get_uniprot_details(unicode, debug):
```

When given a uniprot code, will return a dictionary containing details from 
the uniprot website.

### Parameters

**unicode** : str  
The uniprot code used to fetch uniprot details.

**debug** : bool, optional  
Whether to print notifications.


### Returns

**dict**  
A dictionary relating to different information from uniprot, such as protein names
and descriptions.

### Examples

```python
>>> get_uniprot_details("Q13426")

{'uniprot name': 'DNA repair protein XRCC4',
 'uniprot abbreviation': 'XRCC4_HUMAN',
 'localisations': ['Nucleus', 'Chromosome'],
 'description': 'The XRCC4 gene encodes a novel protein involved in DNA double-strand break repair and V(D)J recombination.',
 'variants': {(0, 'In dbSNP:rs28383138.'): {'original residue': 'S',
   'variation residue': 'C',
   'position': '12'},
  (1,
   'In SSMED; impairs the protein function in DNA double-strand break repair; dbSNP:rs587779351.'): {'original residue': 'W',
   'variation residue': 'R',
   'position': '43'},
  (2, 'In dbSNP:rs28383151.'): {'original residue': 'A',
   'variation residue': 'T',
   'position': '56'},
  (3,
   'In SSMED; impaired ability to repair DNA double-strand breaks.'): {'original residue': 'D', 'variation residue': 'E', 'position': '82'},
  (4, 'In dbSNP:rs28360135.'): {'original residue': 'I',
   'variation residue': 'T',
   'position': '134'},
  (5, 'In dbSNP:rs28360136.'): {'original residue': 'E',
   'variation residue': 'Q',
   'position': '142'},
  (6, 'In SSMED.'): {'original residue': '',
   'variation residue': '',
   'begin': '161',
...
 'domains': {},
 'active sites': {},
 'binding sites': {},
 'mammalian': 1,
 'function': ['DNA ligation involved in DNA repair',
  'double-strand break repair',
  'double-strand break repair via nonhomologous end joining',
  'immunoglobulin V(D)J recombination',
  'positive regulation of ligase activity',
  'positive regulation of phosphatidylserine exposure on apoptotic cell surface',
  'protein localization to site of double-strand break',
  'response to X-ray']}

>>> accessions = get_uniprot_accessions("3ii6")
    for chain in accessions:
        details = get_uniprot_details(accessions[chain])
        print(chain, accessions[chain], details["uniprot name"], details["uniprot abbreviation"])

A Q13426 DNA repair protein XRCC4 XRCC4_HUMAN
B Q13426 DNA repair protein XRCC4 XRCC4_HUMAN
C Q13426 DNA repair protein XRCC4 XRCC4_HUMAN
D Q13426 DNA repair protein XRCC4 XRCC4_HUMAN
X P49917 DNA ligase 4 DNLI4_HUMAN
Y P49917 DNA ligase 4 DNLI4_HUMAN

```
    
## get_equivalentresidue

```python
get_equivalentresidue(resnum : int, seq1 : str, seq2 : str, flanknum : int = 5, placeholder : str = "!", pass_nan : bool = True, debug : bool = False) -> list:
```

Takes the specified residue from sequence 1 and uses alignment to get its number 
in sequence 2.

### Parameters

**resnum** : int  
    The residue number to convert, starts at 1.

**seq1** : str  
    The sequence that residue number is from.

**seq2** : str  
    The sequence to find that residue in.

**flanknum** : int, optional  
    The number of flanking residues to use in the alignment. The default is 5.

**placeholder** : str, optional  
    The placeholder to be used when a residue is missing, such as the beginning 
    and ends of the sequence. Must be one character. The default is "!".

**pass_nan** : bool, optional  
    Whether to pass or raise an exception when given nan in the sequences. Default is True (pass)

**debug** : bool, optional  
    Whether a message should be printed when failing to find a residue. The default is True.

### Returns

**list**   
    0: converted residue number  
    1: alignment score  

### Examples

```python
>>> get_equivalentresidue(2, "ASDF", "FDSA")

[3, 7]

>>> get_equivalentresidue(6, "QWERTYUIOP", "ASDFGHJKLF")

[nan, 0]

>>> get_equivalentresidue(6, "QWERTYUIOP", "ASDFYGHJYYKLQWERTYUIOP")

[18, 11]
```

## get_uniprot_sequence

```python
get_uniprot_sequence(unicode : str, pass_nan : bool = True, pass_no_output = True, debug : bool = True) -> str
```
Get the sequence of a protein from its uniprot accession code.

### Parameters

unicode : str  
    Uniprot accession code.

pass_nan : bool, optional  
    As long as True, input of nan will be passed, otherwise will raise exception. Default is True.

pass_no_output : bool, optional  
    If false, output of "" will raise exception. Default is True.

debug : bool, optional  
    Should this function print as it goes. Default is True.

### Returns

str
    Residue sequence of the protein to which uniprot accession was given.

### Examples

```python

>>> get_uniprot_sequence("Q9UHD2")

'MQSTSNHLWLLSDILGQGATANVFRGRHKKTGDLFAIKVFNNISFLRPVDVQMREFEVLKKLNHKNIVKLFAIEEETTTRHKVLIMEFCPCGSLYTVLEEPSNAYGLPESEFLIVLRDVVGGMNHLRENGIVHRDIKPGNIMRVIGEDGQSVYKLTDFGAARELEDDEQFVSLYGTEEYLHPDMYERAVLRKDHQKKYGATVDLWSIGVTFYHAATGSLPFRPFEGPRRNKEVMYKIITGKPSGAISGVQKAENGPIDWSGDMPVSCSLSRGLQVLLTPVLANILEADQEKCWGFDQFFAETSDILHRMVIHVFSLQQMTAHKIYIHSYNTATIFHELVYKQTKIISSNQELIYEGRRLVLEPGRLAQHFPKTTEENPIFVVSREPLNTIGLIYEKISLPKVHPRYDLDGDASMAKAITGVVCYACRIASTLLLYQELMRKGIRWLIELIKDDYNETVHKKTEVVITLDFCIRNIEKTVKVYEKLMKINLEAAELGEISDIHTKLLRLSSSQGTIETSLQDIDSRLSPGGSLADAWAHQEGTHPKDRNVEKLQVLLNCMTEIYYQFKKDKAERRLAYNEEQIHKFDKQKLYYHATKAMTHFTDECVKKYEAFLNKSEEWIRKMLHLRKQLLSLTNQCFDIEEEVSKYQEYTNELQETLPQKMFTASSGIKHTMTPIYPSSNTLVEMTLGMKKLKEEMEGVVKELAENNHILERFGSLTMDGGLRNVDCL'
```

## reverse_sequence

```python
reverse_sequence(sequence, seq_start = False, seq_end = False)
```

Takes input characters and reverses it. If seq_start and seq_end are given, function
returns where seq_start and seq_end are in the new sequence.

### Parameters

**sequence** : str  
    The sequence to reverse.

**seq_start** : int, optional  
    Region position that will be returned in the new sequence. The default is False.

**seq_end** : int, optional  
    Region position that will be returned in the new sequence. The default is False.

### Returns

**str** or  
    **tuple**  
        [0]: The reversed sequence.  
        [1]: The new start position of the region.  
        [2]: The new end position of the region.

### Examples

    >>> reverse_sequence("QWERTYUIOP")

    'POIUYTREWQ'

    >>> reverse_sequence("QWERTYUIOP", 2, 4)

    ('POIUYTREWQ', 7, 9)

## convert_region

```python
convert_region(start_sequence: str, start_region : Union[int, list], end_sequence : str, debug = False) -> dict
```


Takes a range referencing a region in a start sequence and returns where this region is in an end sequence. Slower than get_equivalentresidue but much more versatile due to ability to convert regions.
Use biological sequence numbers (start at 1)

### Parameters

**start_sequence** : str  
    The sequence which the known region belongs to.

**start_region** : int | list
    The residue or residue range (as a list) that this region occupies.

**end_sequence** : str  
    The new sequence where this region should be detected.

**debug** : bool, optional  
    Should progress be printed.

### Returns

**dict**
    Contains "start" and "end" as the start and end of the sequence, and "score" 
    as the alignment score of the new region.

### Examples

```python
>>> convert_region("QWERTYUIOP", [2,5], "ASDFGHJKQWERTYUIOSDFGHJK")

{'start': 10, 'end': 13, 'score': 75.0}

>>> convert_region("QWERTYUIOP", 7, "ASDFGHJKQWERTYUIOSDFGHJK")

{'start': 15, 'end': 15, 'score': 72.72727272727273}
```

## convert_regions

```python
convert_regions(regions : dict, seq1 : str, seq2) -> dict
```

Takes a dictionary of significant regions as an input, converts them from sequence 1
to sequence 2. If a dict of sequences is given for seq2, then will find the region 
in all those sequences.

### Parameters

**regions** : dict or str in dict format  
    Input set of regions to be converted.

**seq1** : str  
    Sequence the regions are from.

**seq2** : str or dict  
    Sequence the regions should be converted to.

### Returns

**dict**  
    New dictionary of regions that have been converted from seq1 sequence to seq2 
    sequence.

### Examples

```python

# Find where domains from TBK1 appear in similar IKKɛ
>>> tbk1_sequence = get_uniprot_sequence("Q9UHD2")
>>> tbk1_domains = get_uniprot_details("Q9UHD2")["domains"]
>>> ikk_sequence = get_uniprot_sequence("Q14164")
>>> convert_regions(tbk1_domains, tbk1_sequence, ikk_sequence)

{(0, 'Protein kinase'): {'begin': 18, 'end': 305, 'score': 63.005050505050505},
 (1, 'Ubiquitin-like'): {'begin': 310, 'end': 327, 'score': 51.26262626262626}}

```


## get_oximouse_data

```python
get_oximouse_data(age : str)
```
    
 Fetch oximouse data as a dataframe. Age can be "aged", "young", or "detected".
Aged or young will return the data for the aged or young mice, respectively.
Detected returns a list of every cysteine that was detected with oximouse.
Oximouse data maps oxidation of specific cysteines in different body parts in
aged or young mice.

Reference:  
Xiao, H., Jedrychowski, M. P., Schweppe, D. K., Huttlin, E. L., Yu, Q., Heppner, D. E., Li, J., Long, J., Mills, E. L., Szpyt, J., He, Z., Du, G., Garrity, R., Reddy, A., Vaites, L. P., Paulo, J. A., Zhang, T., Gray, N. S., Gygi, S. P., & Chouchani, E. T. (2020). A Quantitative Tissue-Specific Landscape of Protein Redox Regulation during Aging. Cell, 180(5), 968-983.e24. https://doi.org/10.1016/j.cell.2020.02.012


### Parameters

**age** : str  
    "aged", "young", or "detected"

### Returns
Oximouse dataset as a pandas dataframe

### Examples

```python

>>> old = get_oximouse_data("aged")

Downloaded aged oximouse data. Please cite
Xiao, H., Jedrychowski, M. P., Schweppe, D. K., Huttlin, E. L., Yu, Q., Heppner, D. E., Li, J., Long, J., Mills, E. L., Szpyt, J., He, Z., Du, G., Garrity, R., Reddy, A., Vaites, L. P., Paulo, J. A., Zhang, T., Gray, N. S., Gygi, S. P., & Chouchani, E. T. (2020). A Quantitative Tissue-Specific Landscape of Protein Redox Regulation during Aging. Cell, 180(5), 968-983.e24. https://doi.org/10.1016/j.cell.2020.02.012

>>> old

Uniprot ID aged	Gene symbol aged	Site aged	ModScore aged	Motif aged	BAT aged	Brain aged	Epi aged	Heart aged	Kidney aged	...	BAT dev aged	Brain dev aged	Epi dev aged	Heart dev aged	Kidney dev aged	Liver dev aged	Lung dev aged	SKM dev aged	Spleen dev aged	SubQ dev aged
0	A0JNU3	ASPG	503	1000.0	DVGTELCRLASRG	NaN	NaN	NaN	NaN	3.48	...	NaN	NaN	NaN	NaN	0.5	0.21	NaN	NaN	NaN	NaN
1	A0JNU3	ASPG	558	1000.0	SFKDSVCAQPQPH	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	0.16	NaN	NaN	NaN	NaN
2	A0JNU3	ASPG	198	1000.0	RRFAAFCSPNLPP	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	2.66	NaN	NaN	NaN	NaN
3	A0JNU3	ASPG	376	1000.0	QDGMLGCRVAWLL	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	0.21	NaN	NaN	NaN	NaN
4	A1BN54	ACTN1	774	1000.0	TDDFRACLISMGY	0.31	7.02	10.72	7.29	2.7	...	0.29	0.97	0.9	0.87	0.45	0.71	0.43	NaN	9.83	2.55
...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
25924	G5E8Z3	2310050C09RIK	278	0.0	CDHEDDCCCxxxx	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	NaN	0.95	NaN	NaN
25925	P01027	C3	693	0.0	DKGLRKCCEDGMR	NaN	NaN	3.95	NaN	NaN	...	NaN	NaN	3.94	NaN	NaN	NaN	NaN	NaN	NaN	NaN
25926	P10605	CTSB	108	0.0	QGSCGSCWAFGAV	NaN	NaN	NaN	NaN	22.82	...	NaN	NaN	NaN	NaN	1.38	NaN	NaN	NaN	NaN	NaN
25927	Q80W15	IGFBPL1	57	0.0	ARDECGCCARCLG	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	2.03
25928	Q8CG65	SSPO	2307	0.0	SPSQLRCGSGECL	NaN	NaN	NaN	NaN	NaN	...	NaN	NaN	NaN	NaN	NaN	NaN	1.64	NaN	NaN	NaN
25929 rows × 28 columns

```

## get_CIF_structure
```python
get_CIF_structure(pdb_id : str, folder : str = "structures", strict = True, debug = False)
```
Downloads a CIF protein structure from the PDB for a given PDB ID. This should be used 
instead of previous get_PDB_structure()

### Parameters

**pdb_id** : str  
    PDB ID of the structure to download.

**folder** : str, optional  
    Folder to save the structure to. The default is "structures".

**strict** : bool, optional  
    Whether to raise an exception if the structure is not found. The default is True.

**debug** : bool, optional  
    Whether to print messages as it goes. The default is False.
    

### Returns

None. Saves the structure to the specified folder, default "structures".


## get_alphafold_structure

```python
get_alphafold_structure(uniprot_code : str, folder : str = "structures", extension = "ent", strict = False, debug = False)
```

Downloads a structure from AlphaFold for a given uniprot code.

### Parameters

**uniprot_code** : str  
    Uniprot code of the protein to download.

**folder** : str, optional
    Folder to save the structure to. The default is "structures".

**extension** : str, optional  
    File extension to give to the structure. The default is "ent".

**strict** : bool, optional  
    Whether to raise an exception if the structure is not found. The default is False.

**debug** : bool, optional  
    Whether to print messages as it goes. The default is False.

### Returns

None. Saves the AlphaFold structure for that uniprot code.

### Examples 


```python
>>> get_alphafold_structure("Q13426")

Downloading structure for Q13426 from AlphaFold. Please cite: 
Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2
```

## run_propka

```python
run_propka(input_file, structure_folder = "pdb", structure_extension = "ent", propka_folder = "propka/", check = True)
```

Checks if a propka file exists, if not then it attempts to make compute one.

### Parameters

**input_file** : str  
    The name of the file to compute propka for.

**structure_folder** : str, optional  
    The folder to look for the structure in. The default is "pdb".

**structure_extension** : str, optional  
    The extension of the structure file. The default is "ent".

**propka_folder** : str, optional  
    The folder to save the propka file in. The default is "propka/".

**check** : bool, optional  
    Whether to check if the propka file exists before computing it. The default is True.

### Returns

**i** : propka.run.single  
    The propka object. Also saves it to a file.

## check_structure_for_proximal_atoms

```python
check_structure_for_proximal_atoms(structure_file, residue_1, residue_2, atom_1 = "CA", atom_2 = "CA", max_distance = 10)
```

Open a protein structure (or structure) and search for two
residues that are within a specified distance of each other. Returning a list of
dictionaries with two residues and their distance from each other.

### Parameters

**structure_file** : str  
    Path to the structure file

**residue_1** : int  
    Residue number of the first residue

**residue_2** : int  
    Residue number of the second residue

**atom_1** : str, optional  
    Name of the atom in the first residue. Default is "CA"

**atom_2** : str, optional  
    Name of the atom in the second residue. Default is "CA"

**distance** : int, optional  
    Distance cutoff. Default is 10.

### Returns

**list**
    List of dicionaries containing the two residues and their distance from each other.

### Examples

```python
>>> get_alphafold_structure("A2A5R2")
>>> check_structure_for_proximal_atoms("structures/A2A5R2.ent", "CYS", "CYS", atom_1 = "SG", atom_2 = "SG", max_distance = 5)

[{'residue number A': 125,
  'chain A': 'A',
  'residue number B': 76,
  'chain B': 'A',
  'distance': 3.7829864},
 {'residue number A': 514,
  'chain A': 'A',
  'residue number B': 573,
  'chain B': 'A',
  'distance': 3.7547925},
 {'residue number A': 725,
  'chain A': 'A',
  'residue number B': 763,
  'chain B': 'A',
  'distance': 4.4332957},
 {'residue number A': 1350,
  'chain A': 'A',
  'residue number B': 1286,
  'chain B': 'A',
  'distance': 3.520143},
 {'residue number A': 1417,
  'chain A': 'A',
  'residue number B': 1457,
  'chain B': 'A',
  'distance': 4.3076897}]
```

## extract_chain_sequences_from_structure

```python
extract_chain_sequences_from_structure(structure)
```

Extract sequences from each chain in a protein structure.

### Parameters

**structure** : Bio.PDB.Structure  
    The protein structure to extract sequences from.

### Returns

**dict**  
    Dictionary with chain IDs as keys and sequences as values. Gaps and unknown residues are represented with "!".

### Examples

```python
>>> from Bio import PDB
>>> structure = PDB.PDBParser().get_structure("struc", "1abc.pdb")
>>> sequences = extract_chain_sequences_from_structure(structure)
>>> print(sequences)
{'A': 'MKWVTFISLLLLFSSAYS...', 'B': 'VLSPADKTNVKAAW...'}
```

## extract_interactions

```python
extract_interactions(structure, max_distance=4, strict=True, debug=False)
```

Iterate through the chains in a structure and extract regions that interact with ions, ligands, and other chains in the structure.

### Parameters

**structure** : Bio.PDB.Structure  
    Structure to extract interactions from.

**max_distance** : int, optional  
    Maximum distance for an interaction to be considered. The default is 4.

**strict** : bool, optional  
    Whether to raise an exception if no interactions are found. If False then if no interactions are found will return an empty dictionary. The default is True.

**debug** : bool, optional  
    Whether to print debug information. The default is False.

### Returns

**dict**  
    Dictionary containing interaction information for each chain in the structure. The dictionary maps chain IDs to dictionaries of interactions, where each interaction specifies the interacting residue numbers and type of interaction (e.g. ion, ligand, other chain).

### Examples

```python
>>> from Bio import PDB
>>> structure = PDB.PDBParser().get_structure("struc", "3OCP.ent")
>>> extract_interactions(structure)
{'A': {'chain B': [[307, 308], [363, 364]], 'Mg': [[15], [36]]}}
```

## download_structures

```python
download_structures(structures: list, max_concurrent: int = 5, folder: str = "structures", debug: bool = False)
```

Concurrently download multiple protein structures from the PDB.

### Parameters

**structures** : list  
    List of structure IDs (PDB codes) to download.

**max_concurrent** : int, optional  
    Maximum number of concurrent downloads. Default is 5.

**folder** : str, optional  
    Directory where structures will be saved. Default is "structures".

**debug** : bool, optional  
    Whether to print debug information. Default is False.

### Returns

None  
    Structures are downloaded and saved to the specified folder.

### Examples

```python
>>> download_structures(["1ABC", "2DEF", "3GHI"])
# Downloads these three structures concurrently to the "structures" folder
```