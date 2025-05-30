"""
Test-suite for strucseq.
"""

import string
import random
import os
import time
import ast
import numpy as np
from tqdm import tqdm
import strucseq as sq

def mutate_string(s: str) -> str:
    """
    Mutate a string by removing, inserting, or replacing a random letter.
    """
    if not s:
        return s  # avoid crashing on empty string

    idx = random.randrange(len(s))
    rand = random.random()
    if rand < 0.33:
        # Remove the character at idx
        return s[:idx] + s[idx+1:]
    elif rand < 0.67:
        # Insert a random character at idx
        new_char = random.choice(string.ascii_letters)
        return s[:idx] + new_char + s[idx:]
    else:
        # Replace the character at idx with a random character
        new_char = random.choice(string.ascii_letters)
        return s[:idx] + new_char + s[idx+1:]
    
def test_get_equivalentresidue():
    """
    Test the get_equivalentresidue function.
    """
    print("Testing get_equivalentresidue")
    sequence = "MERKISRIHLVSEPSITHFLQVSWEKTLESGFVITLTDGHSAWTGTVSESEISQEADDMAMEKGKYVGELRKALLSGAGPADVYTFNFSKESCYFFFEKNLKDVSFRLGSFNLEKVENPAEVIRELICYCLDTIAENQAKNEHLQKENERLLRDWNDVQGRFEKCVSAKEALETDLYKRFILVLNEKKTKIRSLHNKLLNAAQEREKDIKQEGETAICSEMTADRDPVYDESTDEESENQTDLSGLASAAVSKDDSIISSLDVTDIAPSRKRRQRMQRNLGTEPKMAPQENQLQEKENSRPDSSLPETSKKEHISAENMSLETLRNSSPEDLFDEI"

    # If the result dictionary already exists, read it
    if os.path.exists("test_get_equivalentresidue.txt"):
        with open("test_get_equivalentresidue.txt", "r", encoding="utf-8") as f:
            old_result = ast.literal_eval(f.read())
            old_time = old_result["time"]
            old_result.pop("time")

        # Then test each result
        test_dict = {key : "" for key in old_result}
        
        start_time = time.time()
        for key in tqdm(test_dict):
            test_dict[key] = sq.get_equivalentresidue(165, sequence, key)
        end_time = time.time()
        time_taken = end_time - start_time
        print(f"Time taken: {time_taken} seconds")

        # Check if the results are the same
        if old_result == test_dict:
            print(f"Passed on all tests! Time taken: {time_taken} seconds compared to original {old_time} seconds")
        else:
            print("Test failed")
            for key in test_dict:
                if old_result[key] != test_dict[key]:
                    print(f"Test failed for {key}, got {test_dict[key]} instead of {old_result[key]}")
            raise AssertionError("Test failed")
    else:
        # If the result dictionary does not exist, create it
        print("Building original result dictionary")
        result = {}

        # Make a dictionary of 10000 mutated sequences
        new_sequence = sequence
        test_dict = {sequence : ""}
        for _ in range(10000):
            new_sequence = mutate_string(new_sequence).upper()
            test_dict[new_sequence] = ""

        # Time how long it takes to get the results, and save the results.
        start_time = time.time()
        for key in tqdm(test_dict):
            result[key] = sq.get_equivalentresidue(165, sequence, key)
        end_time = time.time()

        # Print and save the results
        time_taken = end_time - start_time
        print(f"Time taken: {time_taken} seconds")
        result["time"] = time_taken
        with open("test_get_equivalentresidue.txt", "w", encoding="utf-8") as f:
            f.write(str(result))

test_get_equivalentresidue()