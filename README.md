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

## reverse_sequence

```python
reverse_sequence(sequence, seq_start = False, seq_end = False):
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
convert_region(start_sequence: str, start_region : Union[int, list], end_sequence : str, debug = False) -> dict:
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
convert_regions(regions : dict, seq1 : str, seq2) -> dict:
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

# Find if and where modifications from TBK1 that also appear in similar IKKɛ
>>> convert_regions(get_uniprot_details("Q9UHD2")["domains"], "MQSTSNHLWLLSDILGQGATANVFRGRHKKTGDLFAIKVFNNISFLRPVDVQMREFEVLKKLNHKNIVKLFAIEEETTTRHKVLIMEFCPCGSLYTVLEEPSNAYGLPESEFLIVLRDVVGGMNHLRENGIVHRDIKPGNIMRVIGEDGQSVYKLTDFGAARELEDDEQFVSLYGTEEYLHPDMYERAVLRKDHQKKYGATVDLWSIGVTFYHAATGSLPFRPFEGPRRNKEVMYKIITGKPSGAISGVQKAENGPIDWSGDMPVSCSLSRGLQVLLTPVLANILEADQEKCWGFDQFFAETSDILHRMVIHVFSLQQMTAHKIYIHSYNTATIFHELVYKQTKIISSNQELIYEGRRLVLEPGRLAQHFPKTTEENPIFVVSREPLNTIGLIYEKISLPKVHPRYDLDGDASMAKAITGVVCYACRIASTLLLYQELMRKGIRWLIELIKDDYNETVHKKTEVVITLDFCIRNIEKTVKVYEKLMKINLEAAELGEISDIHTKLLRLSSSQGTIETSLQDIDSRLSPGGSLADAWAHQEGTHPKDRNVEKLQVLLNCMTEIYYQFKKDKAERRLAYNEEQIHKFDKQKLYYHATKAMTHFTDECVKKYEAFLNKSEEWIRKMLHLRKQLLSLTNQCFDIEEEVSKYQEYTNELQETLPQKMFTASSGIKHTMTPIYPSSNTLVEMTLGMKKLKEEMEGVVKELAENNHILERFGSLTMDGGLRNVDCL", "MQSTSNHLWLLSDILGQGATANVFRGRHKKTGDLFAIKVFNNISFLRPVDVQMREFEVLKKLNHKNIVKLFAIEEETTTRHKVLIMEFCPCGSLYTVLEEPSNAYGLPESEFLIVLRDVVGGMNHLRENGIVHRDIKPGNIMRVIGEDGQSVYKLTDFGAARELEDDEQFVSLYGTEEYLHPDMYERAVLRKDHQKKYGATVDLWSIGVTFYHAATGSLPFRPFEGPRRNKEVMYKIITGKPSGAISGVQKAENGPIDWSGDMPVSCSLSRGLQVLLTPVLANILEADQEKCWGFDQFFAETSDILHRMVIHVFSLQQMTAHKIYIHSYNTATIFHELVYKQTKIISSNQELIYEGRRLVLEPGRLAQHFPKTTEENPIFVVSREPLNTIGLIYEKISLPKVHPRYDLDGDASMAKAITGVVCYACRIASTLLLYQELMRKGIRWLIELIKDDYNETVHKKTEVVITLDFCIRNIEKTVKVYEKLMKINLEAAELGEISDIHTKLLRLSSSQGTIETSLQDIDSRLSPGGSLADAWAHQEGTHPKDRNVEKLQVLLNCMTEIYYQFKKDKAERRLAYNEEQIHKFDKQKLYYHATKAMTHFTDECVKKYEAFLNKSEEWIRKMLHLRKQLLSLTNQCFDIEEEVSKYQEYTNELQETLPQKMFTASSGIKHTMTPIYPSSNTLVEMTLGMKKLKEEMEGVVKELAENNHILERFGSLTMDGGLRNVDCL")

{(0, 'Protein kinase'): {'begin': 9, 'end': 310, 'score': 100.0},
 (1, 'Ubiquitin-like'): {'begin': 309, 'end': 385, 'score': 100.0}}

```