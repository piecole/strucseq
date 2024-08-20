import strucseq.strucseq as sq
import ast
import numpy as np


print(sq.get_uniprot_accessions("3ii6") == {'A': 'Q13426', 'B': 'Q13426', 'C': 'Q13426', 'D': 'Q13426', 'X': 'P49917', 'Y': 'P49917'})
print(sq.get_uniprot_accessions("6abo") == {'A': 'Q13426', 'B': 'Q0D2I5'})
print(sq.get_uniprot_accessions(np.nan, strict=False))
try:
    print(sq.get_uniprot_accessions(np.nan, strict=True))
    print(False)
except:
    print(True)