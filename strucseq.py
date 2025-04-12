import time
import gzip
from typing import Union
import ast
import os
import glob
import shutil
import csv
import re
import statistics as stats
import requests
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
from Bio.PDB import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB.NeighborSearch import NeighborSearch

try:
    import propka.run as pk
except ModuleNotFoundError:
    print("PROPKA not installed, will not be able to determine pka through PROPKA.")


pdb_search_enabled = False
try:
    from rcsbsearchapi.search import TextQuery as PDBquery
    pdb_search_enabled = True
except ModuleNotFoundError:
    print("rcsbsearchapi not installed, some functions may not work.")

try:
    from tqdm import tqdm
    tqdm.pandas()
except ModuleNotFoundError:
    def tqdm(iterator, *args, **kwargs):
        return iterator
    print("tqdm not installed, progress bars will not work.")

alphabet = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P",
            "Q","R","S","T","U","V","W","X","Y","Z", "DCY"]
threeletter = ["ALA","DCY","CYS","ASP","GLU","PHE","GLY","HIS","ILE","JXX",
               "LYS","LEU","MET","ASN","OXX","PRO","GLN","ARG","SER","THR",
               "MSE","VAL","TRP","TPO","TYR","SEP"]
threetoone = dict(zip(threeletter, alphabet))

# Fix forward slashes in filepaths and create missing folders
def parse_folder(input_folder : str):
    """ Take a folder string and ensure it is formatted correctly. """
    assert isinstance(input_folder,
                      str), "str expected for input_folder, got" + repr(type(input_folder))
    input_folder = input_folder.replace("\\", "/")
    if input_folder[-1] != "/":
        input_folder = input_folder + "/"

    # Create missing folder
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    return input_folder

get_pKa = { #this dictioary defines the pKas of amino acid side chains
    # According to:
    # https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html#prop
    "R": 12.48, "D": 3.65, "C": 8.18, "E": 4.25, "H": 6.00, "K": 10.53,
    "Y": 10.07, "U":1, "X":5.9, "Z":5.6, "DCY":8.18
    # Selenomethionine "U" has been given as pKa 1 as no value is known
    # Phosphoserine and phospohthreonine are from 16018962
    # R-cysteine given same value as L-cysteine
}
get_charge = {
    "R": 1, "D": -1, "E": -1, "H": 1, "K": 1
    }

get_hydrophobic7 = { # All the amino acid hydrophobicities at pH7
    "A":41,
    "C":49,
    "D":-55,
    "E":31,
    "F":97,
    "G":0,
    "H":8,
    "I":99,
    "K":-23,
    "L":97,
    "M":74,
    "N":-28,
    "P":-49, #reference uses pH2 value for some reason
    "Q":-10,
    "R":-14,
    "S":-5,
    "T":13,
    "U":74, # Selenomethonine same as methionine since no value is known
            # for selenomethionine? - check this
    "V":76,
    "W":97,
    "X":0, # Phospohthreonine unknown what value to give
    "Y":63,
    "Z":0, # Phosphothreonine unknown what value to give
    "DCY":49 # R-cysteine same as L cysteine
}

def get_uniprot_accessions(pdb_id : str) -> dict:
    """
    Fetches the UniProt-to-chain mappings for a given PDB ID from the PDBe API.

    Parameters
    ----------
    pdb_id : str
        The 4-character PDB identifier.

    Returns
    -------
    dict: A dictionary containing the mappings of UniProt IDs to chain IDs and residue ranges.
            Returns None if an error occurs.
    """
    # Ensure the PDB ID is in lowercase
    pdb_id = pdb_id.lower()

    # Construct the API URL
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}'

    # Make the HTTP GET request
    response = requests.get(url)
    response.raise_for_status()  # Raises HTTPError for bad responses

    # Parse the JSON response
    data = response.json()

    # Initialize an empty dictionary to store the chain-to-UniProt mapping
    chain_to_uniprot = {}

    # Extract the UniProt mappings
    for pdb_id, pdb_data in data.items():
        uniprot_entries = pdb_data.get('UniProt', {})
        for uniprot_code, uniprot_data in uniprot_entries.items():
            mappings = uniprot_data.get('mappings', [])
            for mapping in mappings:
                chain_id = mapping.get('chain_id')
                if chain_id:
                    chain_to_uniprot[chain_id] = uniprot_code

    return chain_to_uniprot


#   OLD VERSION OF THE FUNCTION:
def iterate_uniprot_accessions_OLD(in_csv : str,
                                   out_csv : str,
                                   delimiter = "\t",
                                   rows : str = 10000000000000,
                                   debug = True):
    """
    Takes an input CSV with a PDBid header and some custom headers to get the chains from
    and returns a csv with the assession codes and unique structures and chains.

    Parameters
    ----------
    in_csv : str
        The the input CSV with PDBids that is used to generate accession codes.
    out_csv : str
        The output CSV used to store fetched accession codes.
    delimiter : str, optional
        Delimiter of the CSV file
    rows : str, optional
        The max index for fetching, usually for testing. The default is 10000000000000.
    debug : bool, optional
        Whether to print progess or other notifications
        
    THIS FUNCTION IS TOO SPECIFIC TO ReDisulphID WRITTING NEW VERSION NEXT

    """

    data = pd.read_csv(in_csv, delimiter = delimiter)
    try:
        previous_file = pd.read_csv(out_csv, delimiter = delimiter)
        fetched_chains = list(previous_file["PDB"].drop_duplicates())
        if debug:
            print("already got:", fetched_chains)
    except:
        previous_file = pd.DataFrame()
        fetched_chains = []
        if debug:
            print("starting new PDB list")

    apd = data[["PDBid", "a chain"]].rename(columns = {"a chain" : "chain"})
    bpd = data[["PDBid", "b chain"]].rename(columns = {"b chain" : "chain"})
    uniquechains = pd.concat([apd,
                              bpd], ignore_index = True).reset_index(drop = True).drop_duplicates()

    uniquePDBs = uniquechains["PDBid"].drop_duplicates().reset_index(drop = True)

    uniprotpd = pd.DataFrame({"PDB" : [], "chain": [], "uniprot": []})
    i = 0
    e = 0

    pbar = tqdm(np.setdiff1d(list(uniquePDBs), fetched_chains))
    for PDBcode in pbar:
        pbar.set_description("Getting uniprot accessions: " + PDBcode)
        if i < rows:

            chaindetails = {}
            tries = 0
            while tries < 3 and not chaindetails: # If no details are returned,
                                                  # try a couple more times
                # This is because sometimes no details are returned when there
                # should be some, probably due to anti-scraping measures
                chaindetails = get_uniprot_accessions(PDBcode)
                tries = tries + 1

            if debug:
                if tries == 3:
                    print("no accessions found")
                if 1 < tries < 3:
                    print("found accessions in", str(tries), "tries.")
            if debug:
                print(chaindetails)

            for letter in chaindetails: # Iterate each letter in the dictionary
                #   Create a new row with the ID, chain letter, and uniprot code
                uniprotpd.loc[e] = [PDBcode, letter, chaindetails[letter]]
                e = e + 1
        i = i + 1

    uniprotpd = pd.concat([uniprotpd, previous_file])
    try:
        uniprotpd.to_csv(out_csv, sep=delimiter, index = False)
    except:
        if debug:
            print("Failed save, using utf-8 instead.")
        uniprotpd.to_csv(out_csv, sep=delimiter, encoding='utf-8', index = False)

def iterate_uniprot_accessions(in_csv : str,
                               chain_cols,
                               out_csv : str,
                               delimiter = "\t",
                               debug = True):
    """
    Takes an input CSV with a PDBid header and specified custom header(s) to get the
    uniprot ID of the chains and saves a CSV with the assession codes and unique structures
    and chains.

    Parameters
    ----------
    in_csv : str
        The the input CSV with PDBids that is used to generate accession codes.
    chain_cols : list or str
        Column or list of columns that should be used to identify the chain(s) when
        fetching their accessions.
        e.g. ["a chain", "b chain"]
    out_csv : str
        The output CSV used to store fetched accession codes.
    delimiter : str, optional
        Delimiter of the CSV file
    debug : bool, optional
        Whether to print progess or other notifications

    """

    data = pd.read_csv(in_csv, delimiter = delimiter)
    try:
        previous_file = pd.read_csv(out_csv, delimiter = delimiter)
        fetched_chains = list(previous_file["PDB"].drop_duplicates())
        if debug:
            print("already got:", fetched_chains)
    except:
        previous_file = pd.DataFrame()
        fetched_chains = []
        if debug:
            print("starting new PDB list")

    #   Make separate dataframes for each of the chain cols and then combine them
    uniquechains = pd.DataFrame()

    #   Parse chain_cols into a list if it is a string
    if isinstance(chain_cols, str):
        chain_cols = [chain_cols]

    for chain_col in chain_cols:
        #   Make the df
        new_df = data[["PDBid", chain_col]].rename(columns = {chain_col : "chain"})
        #   Combine it with the previously made df
        uniquechains = pd.concat([uniquechains,
                                  new_df], ignore_index = True).reset_index(drop = True).drop_duplicates()

    uniquePDBs = uniquechains["PDBid"].drop_duplicates().reset_index(drop = True)
    
    uniprotpd = pd.DataFrame({"PDB" : [], "chain": [], "uniprot": []})
    e = 0

    pbar = tqdm(np.setdiff1d(list(uniquePDBs), fetched_chains))
    for PDBcode in pbar:
        pbar.set_description("Getting uniprot accessions: " + PDBcode)

        chaindetails = {}
        tries = 0
        while tries < 3 and not chaindetails: # If no details are returned,
                                              # try a couple more times
            # This is because sometimes no details are returned when there
            # should be some, probably due to anti-scraping measures
            chaindetails = get_uniprot_accessions(PDBcode)
            tries = tries + 1

        if debug:
            if tries == 3:
                print("no accessions found")
            if 1 < tries < 3:
                print("found accessions in", str(tries), "tries.")
        if debug:
            print(chaindetails)

        for letter in chaindetails: #   Iterate each letter in the dictionary
            #   Create a new row with the ID, chain letter, and uniprot code
            uniprotpd.loc[e] = [PDBcode, letter, chaindetails[letter]] 
            e = e + 1


    uniprotpd = pd.concat([uniprotpd, previous_file])
    try:
        uniprotpd.to_csv(out_csv, sep=delimiter, index = False)
    except:
        if debug:
            print("Failed save, using utf-8 instead.")
        uniprotpd.to_csv(out_csv, sep=delimiter, encoding='utf-8', index = False)

def get_uniprot_details(unicode : str, debug = False) -> dict:
    """
    
    When given a uniprot code, will return a dictionary containing details from 
    the uniprot website.

    Parameters
    ----------
    unicode : str
        The uniprot code used to fetch uniprot details.
    debug : bool, optional
        Whether to print notifications.

    Returns
    -------
    dict
        A dictionary relating to different information from uniprot, such as protein names
        and descriptions.

    """

    try:
        url = "https://www.uniprot.org/uniprot/" + unicode + ".xml"
        if debug:
            print("downloading", url)
    except:
        raise Exception("Give a uniprot accession code, recieved " + repr(unicode))

    output = {}
    if debug:
        print(f"Extracting Uniprot information from code: {unicode}.")

    # Fetch the xml version of the uniprot site in beautifulsoup
    # Added try loop due to error: ConnectionError: ('Connection aborted.',
    # RemoteDisconnected('Remote end closed connection without response'))
    tries = 0
    while tries < 10:
        try:
            if debug:
                print(f"Accessing {url}")
            soup = BeautifulSoup(requests.get(url).text, "lxml")
            # If successful break out of while loop
            break
        except (ConnectionResetError,
                requests.exceptions.ProtocolError,
                requests.exceptions.ConnectionError,
                requests.exceptions.SSLError):
            tries += 1
            print(f"Connection error while fetching uniprot data, waiting {tries**2} seconds before retrying...")
            time.sleep(tries**2)
    else:
        print("Failed to fetch uniprot data repeatedly.")
        soup = BeautifulSoup(requests.get(url).text, "lxml")

    
    #   Get the protein name.
    try:
        output["uniprot name"] = soup.find_all("fullname")[0].text
    except IndexError:
        output["uniprot name"] = np.nan
    
    #   Get the uniprot abbreviation.
    try:
        output["uniprot abbreviation"] = soup.find_all("name")[0].text
    except IndexError:
        output["uniprot abbreviation"] = np.nan
    
    try:
        #   Search for the list of localisations and then put each one in a list
        localisations = [x.text for x in soup.find("comment",
                                                   {"type" : "subcellular location"}).find_all("location")]
        #   Output this list
        output["localisations"] = localisations
    except:
        #   If there is no such localisation
        output["localisations"] = np.nan
    
    #   Get the first description
    try:
        output["description"] = soup.find("title").text
    except:
        #   If no description found.
        output["description"] = np.nan
        
    #   Get disease variants
    variants = {}
    for variant in soup.find_all("feature", {"type":"sequence variant"}):
        try:
            variant["description"] # Check the variant has a description
        except KeyError:
            variant["description"] = "none"
        try: #  Check whether deletions
            variant.find("original").text
            variant.find("variation").text
        except: #   If they are deletions (shown by no residue change)
            try:
                variant.append(soup.new_tag("original", text = variant.find("position")["position"])) # Make the original residue that res nubmer
            except: #   If its muliple residues then make that original residue a resnumber range
                variant.append(soup.new_tag("original", text = variant.find("begin")["position"] + "-" + variant.find("end")["position"]))
            variant.append(soup.new_tag("variation", text = "deletion")) #make the new residue a deletion
        try:
            #   Single point change
            variants[len(variants), variant["description"]] = {"original residue" : variant.find("original").text, "variation residue" : variant.find("variation").text, "position" : variant.find("location").find("position")["position"]}
        except:
            #   Multi point changes
            variants[len(variants), variant["description"]] = {"original residue" : variant.find("original").text, "variation residue" : variant.find("variation").text, "begin" : variant.find("location").find("begin")["position"], "end" : variant.find("location").find("end")["position"]}
    output["variants"] = variants
    
    #   Get functional mutations (features)
    mutations = {}
    for mutation in soup.find_all("feature", {"type":"mutagenesis site"}):
        try: #  Check whether deletions
            mutation.find("original").text
            mutation.find("variation").text
        except:
            try: #  If they are deletions (shown by no residue change)
                mutation.append(soup.new_tag("original", text = mutation.find("position")["position"]))
            except: #   If its muliple residues then make that original residue a resnumber range
                mutation.append(soup.new_tag("original", text = mutation.find("begin")["position"] + "-" + mutation.find("end")["position"]))
            mutation.append(soup.new_tag("variation", text = "deletion")) #make the new residue a deletion
        try:
            #   Single point mutations
            mutations[len(mutations), mutation["description"]] = {"original residue" : mutation.find("original").text, "variation residue" : mutation.find("variation").text, "position" : mutation.find("location").find("position")["position"]}
        except:
            #   Sequence variants
            mutations[len(mutations), mutation["description"]] = {"begin" : mutation.find("location").find("begin")["position"], "end" : mutation.find("location").find("end")["position"]}
    output["functional mutations"] = mutations
    
    #   Get modified residues
    modifications = {}
    for modification in soup.find_all("feature", {"type" : "modified residue"}):
        modifications[len(modifications), modification["description"]] = {"position" : modification.find("location").find("position")["position"]}
    output["modifications"] = modifications
            
    #   Get functional sites/region of interest
    regions = {}
    for region in soup.find_all("feature", {"type":["region of interest", "short sequence motif"]}):
        try:
            #   Single residue
            regions[len(regions),region["description"]] = {"position" : region.find("location").find("position")["position"]}
        except:
            #   Multiple residue
            for location in ["begin", "end"]:
                try: #  Test if the begining and end have positions
                    _ = region.find(location)["position"]
                except: #   Otherwise assign a new position from the status (likely "unknown")
                    region.find(location)["position"] = region.find(location)["status"]
            regions[len(regions),region["description"]] = {"begin" : region.find("location").find("begin")["position"], "end" : region.find("location").find("end")["position"]}
    output["regions"] = regions
        
    #   Get family & domains
    domains = {}
    for domain in soup.find_all("feature", {"type":"domain"}):
        try:
            #   Single residue domain
            domains[len(domains), domain["description"]] = {"position" : domain.find("location").find("position")["position"]}
        except:
            #   Multi residue domain
            domains[len(domains), domain["description"]] = {"begin" : domain.find("location").find("begin")["position"], "end" : domain.find("location").find("end")["position"]}
    output["domains"] = domains
    
    #   Active site
    active_sites = {}
    for active_site in soup.find_all("feature", {"type":["active site"]}):
        try:
            active_site["description"] #    Test whether the active site has a description
        except:
            active_site["description"] = "none"

        active_sites[len(active_sites),active_site["description"]] = {"position" : active_site.find("location").find("position")["position"]}
            
    output["active sites"] = active_sites
    
    #   Binding site
    binding_sites = {}
    for binding_site in soup.find_all("feature", {"type" : ["binding site"]}):
        try:
            #   Single binding residue
            binding_sites[len(binding_sites),binding_site.find("ligand").find("name").text] = {"position" : binding_site.find("location").find("position")["position"], "reference" : binding_site.find("ligand").find("dbreference")["id"]}
            try:
                binding_sites[len(binding_sites),binding_site.find("ligand").find("name").text] = []
                binding_sites[len(binding_sites),binding_site.find("ligand").find("name").text] = {"position" : binding_site.find("location").find("position")["position"]} #if there is no DB reference
            except:
                #   Binding sequence
                binding_sites[len(binding_sites),binding_site.find("ligand").find("name").text] = {"begin" : binding_site.find("location").find("begin")["position"], "end" : binding_site.find("location").find("end")["position"]}
        except:
            pass
            
    output["binding sites"] = binding_sites
    
    
    #   Default mammalian to 0
    output["mammalian"] = 0
    #   Find all taxons and convert them to a list, then check if one of them is mammalian
    if "Mammalia"  in [x.text for x in soup.find_all("taxon")]:
        output["mammalian"] = 1
    
    #   Get all the terms where processes are mixed up with component and molecular function
    processes = [x.get("value") for x in soup.find_all("property", {"type" : "term"})]
    #   Strip out only the processes
    output["function"] = [x.split(":")[1] for x in processes if "P:" in x]
     
    if debug:
        for i in output:
            print(i, output[i])
    return output

#OLD TOO SPECIFIC CODE:
def iterate_uniprot_details_OLD(in_csv : str, out_csv : str, uniprot_csv : str = None, species : str = None, debug = True):
    """
    Takes an input CSV, uses chain uniprot accessions to add uniprot data to a screen
    product, which is then output as another CSV.

    Parameters
    ----------
    in_csv : str
        The filepath to an input CSV, which must contain columns "uniprot a" and 
        "uniprot b".
    out_csv : str
        The intended filepath of the output CSV.
    uniprot_csv : str, optional
        If a uniprot CSV is then specified use that. The default is None, in which 
        case a column in the input CSV can be used. NEEDS CONFIRMATION
    species : str, optional
        Can specify a species if there is a species-specific column of uniprots 
        that should be used for fetching data. The default is None.
        If a species is specified then it will look in that column for the uniprot 
        accession, rather than fetching the uniprot accession from the separate CSV.
    debug : str, optional
        Whether to print progess.

    Returns
    -------
    A CSV file of the input database populated with protein information from uniprot.
    
    New columns will be the species, then "a chain", or "b chain" followed by:
        'uniprot name', 'uniprot abbreviation', 'localisations', 'description', 'variants',
        'functional mutations', 'modifications', 'regions', 'domains', 'active sites', 
        'binding sites', 'mammalian', and 'function'.
        
    THIS FUNCTION IS TOO SPECIFIC TO ReDisulphID

    """
    
    #   Converts the screen product and uniprot csvs to pandas dataframes
    data = pd.read_csv(in_csv, delimiter = "\t")
    if uniprot_csv != None:
        uniprotdata = pd.read_csv(uniprot_csv, delimiter = "\t")
    
    chains = ["a", "b"]
    #   Combines the known uniprot accessions into the screen product dataframe
    for chain in chains: #  Adds the uniprot based on on 'a chain' then 'b chain'
        if uniprot_csv is not None: #   If uniprotdata was specified, then add this to the dataframe
            data = pd.merge(data, uniprotdata, how = "left", left_on=["PDBid", chain + " chain"], right_on=["PDB", "chain"]) #  Looks at the chain letter and pdb accession code and adds the row from the uniprot dataframe
            data = data.drop(["PDB", "chain"], axis = 1) #  Removes the added PDB code and chain since this information is duplicate
            if species is None:
                #   Renames the new 'chain' column to specify a or b
                data = data.rename(columns = {"uniprot" : "uniprot " + chain}) #    
            else:
                data = data.rename(columns = {"uniprot" : chain + " " + species + " accession code"})
    
    data = data.loc[:,~data.columns.duplicated()].copy()
    
    #if a species is specified, then look in that column for the accession code otherwise look in default column
    accession_columns = ["a", "b"]
    if species is None:
        accession_columns = ["uniprot " + x for x in chains] #default column
    else:
        try:
            accession_columns = [x + " " + species + " accession code" for x in chains]
        except:
            raise Exception("Expected species as string, got " + repr(species))
    
    #accesses the internet using the new uniprots to get uniprot data using get_uniprot_details()
    #creates a new series of unique uniprot codes

    uniqueuniprots = pd.concat([data["uniprot a"], data["uniprot b"]],
                               ignore_index = True).reset_index(drop = True).drop_duplicates()
    uniqueuniprots = pd.DataFrame(uniqueuniprots) #turns this series into a 1 column dataframe
    uniqueuniprots.columns = (["uniprot"]) #gives the 1 column the name 'uniprots'
    uniqueuniprots = uniqueuniprots.dropna() #drops any empty enties
    #drops entries that say "no accession code" (which were returned before when they couldn't be found)
    uniqueuniprots = uniqueuniprots[uniqueuniprots["uniprot"] != "no accession code"]
    
    #creates a test data row based on XRCC4 so we know what the columns are called   
    testset = get_uniprot_details("Q13426")
    testset["uniprot"] = "Q13426"
    dataset = pd.DataFrame(columns = testset.keys()) #creates a new empty dataframe with those columns
    
    #iterates through the uniprots and gets the extra information from the internet added to the new dataset dataframe
    dex = 0
    for index, row in tqdm(uniqueuniprots.iterrows(), total = uniqueuniprots.shape[0]): #iterates through the unique uniprots
        #print(dex, "/", uniqueuniprots.size) #indicates the progress
        dex = dex + 1
        
        if len(row["uniprot"]) > 4: # Check the uniprot code is valid
            # Get uniprot information as a dictionary
            newinfo = get_uniprot_details(row["uniprot"])
            # Add the uniprot accession code to that dictionary
            newinfo["uniprot"] = row["uniprot"]
            # Create a 1 row dataframe of that dictionary, with the keys as the columns 
            for i in newinfo:
                print(f"{i}: {newinfo[i]}")
            newinfo = pd.DataFrame([newinfo],  columns = newinfo.keys())
            # Combine that with the previously constructed dataframe
            dataset = pd.concat([dataset, newinfo], ignore_index = True) # adds 1 row to the dataset dataframe
        else:
            if debug:
                print("not a uniprot accession:", row["uniprot"])
            
            
    #convert species column into dictionary for better access
    accession_columns = {"a chain": accession_columns[0], "b chain" : accession_columns[1]}
    if species is None: #convert the species into text that can be used to name columns
        species = ""
    else:
        species = species + " "
    
    #from the uniprots given by each chain, adds the information from the unique uniprot list to the screen data

    for chain in ["a chain", "b chain"]: #iterate between the two chains
        data = pd.merge(data, dataset, how="left", left_on=[accession_columns[chain]], right_on=["uniprot"]).reset_index(drop = True) #merges the uniprot data to the screen based on the chain uniprots
        data = data.drop(["uniprot"], axis = 1) #drop the new uniprot column since it is a duplicate
        data = data.rename(columns = dict(zip(newinfo.keys(), species + chain + " " + newinfo.keys()))) #renames the new columns to specify their chains
    
    #save all this as a new CSV 
    data.to_csv(out_csv, sep="\t", index = False)
    
def iterate_uniprot_details(in_csv : str,
                            chain_cols : list,
                            out_csv : str,
                            delimiter : str = "\t",
                            uniprot_csv : str = None,
                            species : str = None,
                            debug = True):
    """
    Takes an input pandas DataFrame or CSV file, uses uniprot accessions and provided chain column(s)
    to add uniprot information columns to the data, which is then output as another CSV.

    Parameters
    ----------
    in_csv : str
        The filepath to an input CSV, which must contain columns "uniprot a" and 
        "uniprot b".
    chain_cols : list or str
        String or list of strings specifying which column to get the identity of the 
        peptide chain from.
        E.g. ["a chain", "b chain"]
    out_csv : str
        The intended filepath of the output CSV.
    delimiter : str, optional
        Delimiter that is used in the input CSV, and should be used in the output CSV.
    uniprot_csv : str, optional
        If a uniprot CSV is then specified use that. The default is None, in which 
        case a column in the input CSV can be used. NEEDS CONFIRMATION
    species : str, optional
        Can specify a species if there is a species-specific column of uniprots 
        that should be used for fetching data. The default is None.
        If a species is specified then it will look in that column for the uniprot 
        accession, rather than fetching the uniprot accession from the separate CSV.
    debug : str, optional
        Whether to print progess.

    Returns
    -------
    A CSV file of the input database populated with protein information from uniprot.
    
    New columns will be the species, then specified chain columns followed by:
        'uniprot name', 'uniprot abbreviation', 'localisations', 'description', 'variants',
        'functional mutations', 'modifications', 'regions', 'domains', 'active sites', 
        'binding sites', 'mammalian', and 'function'.
        
    EDITING TO BE NON SPECIFIC
    NEED TO TEST AND FIX THIS WITH pd.DataFrame input

    """

    #   Converts the screen product and uniprot csvs to pandas dataframes
    if isinstance(in_csv, str):
        data = pd.read_csv(in_csv, delimiter = delimiter)
    if isinstance(uniprot_csv, str):
        uniprot_csv = pd.read_csv(uniprot_csv, delimiter = delimiter)
        
    #   Parse chain_cols
    if isinstance(chain_cols, str):
        chain_cols = [chain_cols]
    
    #   Combines the known uniprot accessions into the screen product dataframe
    for chain_col in chain_cols: #  Adds the uniprot based on on 'a chain' then 'b chain'
        if uniprot_csv is not None: #   If uniprot_csv was specified, then add this to the dataframe
            data = pd.merge(data, uniprot_csv, how = "left", left_on=["PDBid", chain_col], right_on=["PDB", "chain"]) #  Looks at the chain letter and pdb accession code and adds the row from the uniprot dataframe
            data = data.drop(["PDB", "chain"], axis = 1) #  Removes the added PDB code and chain since this information is duplicate
            if species is None:
                #Renames the new 'chain' column to specify the chain
                data = data.rename(columns = {"uniprot" : "uniprot " + chain_col}) #    
            else:
                data = data.rename(columns = {"uniprot" : chain_col + " " + species + " accession code"})
    
    # Copy to a new DataFrame that includes only the non-duplicated
    # columns from the original DataFrame.
    data = data.loc[:,~data.columns.duplicated()].copy()
    
    # If a species is specified, then look in that column for the accession
    # code, otherwise look in default column defined by the species.
    accession_columns = ["a", "b"]
    if species is None:
        accession_columns = ["uniprot " + x for x in accession_columns] #default column
    else:
        assert isinstance(species, str), "Expected species as string, found " + repr(species)
        accession_columns = [x + " " + species + " accession code" for x in accession_columns]
    
    #accesses the internet using the new uniprots to get uniprot data using get_uniprot_details()
    #creates a new series of unique uniprot codes
            
    uniqueuniprots = pd.concat([data["uniprot a"], data["uniprot b"]],
                               ignore_index = True).reset_index(drop = True).drop_duplicates()
    uniqueuniprots = pd.DataFrame(uniqueuniprots) #turns this series into a 1 column dataframe
    uniqueuniprots.columns = (["uniprot"]) #gives the 1 column the name 'uniprots'
    uniqueuniprots = uniqueuniprots.dropna() #drops any empty enties
    #drops entries that say "no accession code" (which were returned before when they couldn't be found)
    uniqueuniprots = uniqueuniprots[uniqueuniprots["uniprot"] != "no accession code"]
    
    #creates a test data row based on XRCC4 so we know what the columns are called   
    testset = get_uniprot_details("Q13426")
    testset["uniprot"] = "Q13426"
    dataset = pd.DataFrame(columns = testset.keys()) #creates a new empty dataframe with those columns
    
    # Iterates through the uniprots and gets the extra information from
    # the internet added to the new dataset dataframe
    dex = 0
    for index, row in tqdm(uniqueuniprots.iterrows(),
                           total = uniqueuniprots.shape[0],
                           desc = "Iterating through uniprot site for uniprot details."): #iterates through the unique uniprots
        #print(dex, "/", uniqueuniprots.size) #indicates the progress
        dex = dex + 1
        
        if len(row["uniprot"]) > 4:
            newinfo = get_uniprot_details(row["uniprot"]) #gets uniprot information as a dictionary
            newinfo["uniprot"] = row["uniprot"] #adds the uniprot accession code to that dictionary
            #creates a 1 row dataframe of that dictionary, with the keys as the columns 
            newinfo = pd.DataFrame([newinfo],  columns = newinfo.keys())
            dataset = pd.concat([dataset, newinfo], ignore_index = True) # adds 1 row to the dataset dataframe
        else:
            if debug:
                print("not a uniprot accession:", row["uniprot"])
            
            
    #convert species column into dictionary for better access
    accession_columns = {"a chain": accession_columns[0], "b chain" : accession_columns[1]}
    if species is None: #convert the species into text that can be used to name columns
        species = ""
    else:
        species = species + " "
    
    #from the uniprots given by each chain, adds the information from the unique uniprot list to the screen data

    for chain in ["a chain", "b chain"]: #iterate between the two chains
        data = pd.merge(data, dataset, how="left", left_on=[accession_columns[chain]], right_on=["uniprot"]).reset_index(drop = True) #merges the uniprot data to the screen based on the chain uniprots
        data = data.drop(["uniprot"], axis = 1) #drop the new uniprot column since it is a duplicate
        data = data.rename(columns = dict(zip(newinfo.keys(), species + chain + " " + newinfo.keys()))) #renames the new columns to specify their chains
    
    #save all this as a new CSV
    data.to_csv(out_csv, sep="\t", index = False)

def extract_chain_sequences_from_structure(structure : Structure):
    chain_sequences = {}
    for chain in structure[0]:
        length = max([i.id[1] for i in chain]) #determine how long the chain actually is, ignoring gaps
        if length > 0:
            res_list = ["!" for _ in range(length)] #creating reslist by starting with blank ! marks
            for index, residue in enumerate(chain):
                residue.newresnum = index
                
                #populating res_list, which has gaps
                resname = residue.get_resname()
                try:
                    if resname != "HOH" and residue.id[1] >= 0: #check residue is not water and has a seq number of 0 or more
                        res_list[residue.id[1] -1] = threetoone[resname] #compile chain sequence
                except KeyError:
                    try:
                        res_list[residue.id[1] -1] = "!" #otherwise add exclamation marks
                    except IndexError as e:
                        print("res_list", res_list)
                        print("residue", residue)
                        print("length", length)
                        raise e
                except IndexError:
                    print("res_list:", res_list, "residue:", residue, "resname:", resname)
            chain_sequences[chain.id] = "".join(res_list) #join the res_list compiled previously into strings add to a dictionary with the chain letters as keys
    return chain_sequences

def get_flanking_info(PDB_file : str,
                      amino_acid : str,
                      debug : bool = False) -> tuple:
    """
        
    Takes a PDB file and returns flanking information for all the cysteines and also returns the real sequences of each chain.

    Parameters
    ----------
    PDB_file : str
        The path to the input PDB file.
    amino_acid : str
        The amino acid to list flanking info for, e.g. "CYS".
    debug : bool, optional
        Should this function print its output each time.

    Returns
    -------
    tuple
        [0]: Dataframe with information about flanking residues including sequence, pKas, hydrophobicities, 
        and charges at pH 7.
        
        [1]: Sequences of the chains in the structure.
        
    """
    #catch input errors
    assert isinstance(PDB_file, str), "str expected for PDB_file, found '" + repr(PDB_file) + "' which is " + repr(type(PDB_file))
    assert isinstance(debug, bool), "bool expected for debug, found " + repr(type(debug))
    
    cysteine_list = []
    with gzip.open(PDB_file.encode("unicode_escape"), "rt") as unzipped: #open the structure
        try:
            structure = PDBParser(QUIET=debug).get_structure("struc", unzipped) #parse the structure
        except OSError:
            raise Exception(f"Failed to parse structure from '{PDB_file}'.")
        
        chain_sequences = {}
        for chain in structure[0]: #iterate through chains
            realreslist = [] #make a list to store the residues in a chain
            length = max([i.id[1] for i in chain]) #determine how long the chain actually is, ignoring gaps
            res_list = ["!" for i in range(length)] #creating reslist by starting with blank ! marks
            for index, residue in enumerate(chain):
                #populating realreslist
                residue.newresnum = index
                realreslist.append(residue)
                
                #populating res_list, which has gaps
                resname = residue.get_resname()
                try: 
                    if resname != "HOH" and residue.id[1] >= 0: #check residue is not water and has a seq number of 0 or more
                        res_list[residue.id[1] -1] = threetoone[resname] #compile chain sequence
                except:
                    if debug:
                        print("res list:", res_list, "residue:", residue)
                    res_list[residue.id[1] -1] = "!" #otherwise add exclamation marks
            for residue in realreslist:
                if residue.get_resname() == amino_acid:
                    flanks = []
                    for offset in range(-5, 6):
                        try:
                            flanks.append(threetoone[realreslist[residue.newresnum + offset].get_resname()])
                        except:
                            flanks.append("!")
                    pKas = []
                    hyd = []
                    charge = []
                    for flank in flanks:
                        try:
                            pKas.append(get_pKa[flank])
                        except:
                            pKas.append("")
                        try:
                            hyd.append(get_hydrophobic7[flank])
                        except:
                            hyd.append("")
                        try:
                            charge.append(get_charge[flank])
                        except:
                            charge.append("")
                    cysteine_list.append({"PDBid" : PDB_file.split("pdb")[-1].split(".ent")[0],
                                          "chain" : chain.id,
                                          "residue" : residue.id[1],
                                          "flanking residues" : flanks,
                                          "flanking pKas" : pKas,
                                          "flanking hydrophobicities": hyd,
                                          "flanking charges at pH 7": charge})
                    
            chain_sequences[chain.id] = "".join(res_list) #join the res_list compiled previously into strings add to a dictionary with the chain letters as keys
    
    newdata = pd.DataFrame() #new empty dataframe to start building
    for cysteine in cysteine_list:
        pdrow = pd.DataFrame([[cysteine["PDBid"], cysteine["chain"], cysteine["residue"]]], columns = ["PDBid","chain", "residue"])
        addlettersto = []
        for section in ["flanking residues", "flanking pKas", "flanking hydrophobicities", "flanking charges at pH 7"]: #iterates through the lists in the dictionary created by get_flanking_info()
                for i in range(len(cysteine[section])): #iterates through each list
                    #print("i:", i, "- section:", section, "- newrow:", newrow)
                    pdrow[str(i - 5) + " " + section] = cysteine[section][i] #adds each bit of information on the list as its own column in the empty dataframe, so there is one row
                    addlettersto.append(str(i - 5) + " " + section)
        newdata = pd.concat([newdata, pdrow], ignore_index = True).reset_index(drop = True) #concats this new row into a dataframe to build it
    #print(chain_sequences)
    return newdata, chain_sequences #output the flanking info as a new dataframe for each cysteine, output the chain sequences as a dicitonary

def get_residues(residue : str,
                 flanknum : int,
                 sequence : str,
                 placeholder : str = "!",
                 frameshift : bool = False,
                 strict = False) -> dict:
    """
    Takes a sequence. Creates a dictionary with the residue number and flanking residues 
    of an amino acid that can be found in the sequence.

    Parameters
    ----------
    residue : str
        The residue to return a list of. E.g. "C".
    flanknum : int
        Number of residues to return in the sequences wither side of the desired residue.
    sequence : str
        Sequence in which to find the residues.
    placeholder : str, optional
        Placeholder where there is empty space, such as after the end of a sequence. 
        The default is "!".
    frameshift : bool, optional
        Whether to include also returning frameshifted versions of the sequence, to 
        make sequence alignment more robust. The default is False.
        Frameshift is a bit of a misnomer, just mean (for example) if a residue has 
        been swapped with another next to it, like ABCDEFG compared to ABCEDFG which 
        would otherwise lead to not detecting that residue at all.
    strict : bool, optional
        Whether to raise an exception if the input sequence isn't valid.

    Returns
    -------
    dict
        A dictionary containing sequence numbers as keys with the data attached to each 
        key being a sequence around that residue.
        If frameshift is True, some entries will be frameshifted sequences with a key
        as the residue number followed by +/- and a number denoting the number of 
        frameshifts.

    """

    if strict == True:
        assert isinstance(sequence, str), "Expected str for sequence, got '" + repr(sequence) + "' which is " + repr(type(sequence))

    # The residue number is the actual number rather the list-index number
    reslist = {} # Make the dictionary to fill up
    sequence = "".join([placeholder for i in range(flanknum)]) + sequence

    for position, letter in enumerate(str(sequence)): #iterate through the letters in the sequence
        if(letter == residue): #check that the letter is a cysteine? allowing return of any residue caused problems
            flank = []
            for offset in range(0 - flanknum, flanknum + 1): #count through the residues around the residue
                try:
                    flank.append(sequence[position + offset]) #try to add those flanking residues to 'flank'
                except:
                    flank.append(placeholder) #what to add when there isn't a residue
            reslist[position - flanknum + 1] =  "".join(flank) #add the residue and flanks to the list
    
    #   Add new potential frameshifted residues        
    if frameshift != False:
        shifted = {}
        if isinstance(frameshift, int):
                frameshift += 1
        for res in reslist:    #   For every seq1 residue
        
            #   Set the frameshift level to the flanknum if frameshift not specified
            if isinstance(frameshift, bool):
                frameshift = flanknum + 1

            #   Frameshift every potential residue on the left side
            #   For each number of frameshifts up to the flanknum
            for shift in range(1, frameshift):
                #   Make a new entry to store the shifted sequence as a list
                shifted[str(res) + "-" + str(shift)] = list(reslist[res])
                #   Add an X for every frameshift and remove a residue to compensate
                for i in range(shift):
                    shifted[str(res) + "-" + str(shift)].insert(flanknum, "X")
                    shifted[str(res) + "-" + str(shift)].pop(0)
                    
            #   Frameshift every potential residue on the right side
            for shift in range(1, frameshift):
                shifted[str(res) + "+" + str(shift)] = list(reslist[res])
                for i in range(shift):
                    shifted[str(res) + "+" + str(shift)].insert(flanknum + 1, "X")
                    shifted[str(res) + "+" + str(shift)].pop()
        
        #   Concatonate the shifted sequences into strings so they are the same as 
        #   normal ones.
        for res in shifted:
            shifted[res] = "".join(shifted[res])
        #   Add the shifted sequences to the normal ones
        reslist.update(shifted)
        
    return reslist

amino_acid_groups = {False: [], "positive" : ["R", "H", "K"], "negative": ["D", "E"], "polar": ["S", "T", "N", "Q"], "hydrophobic":["A", "V", "I", "L", "M", "F", "W"]}
def get_equivalentresidue(resnum : int,
                          seq1 : str,
                          seq2 : str,
                          flanknum : int = 5,
                          placeholder : str = "!",
                          pass_nan : bool = True,
                          debug : bool = False) -> list:
    """
    Takes the specified residue from sequence 1 and uses alignment to get its number 
    in sequence 2.

    Parameters
    ----------
    resnum : int
        The residue number to convert, starts at 1.
    seq1 : str
        The sequence that residue number is from.
    seq2 : str
        The sequence to find that residue in.
    flanknum : int, optional
        The number of flanking residues to use in the alignment. The default is 5.
    placeholder : str, optional
        The placeholder to be used when a residue is missing, such as the beginning 
        and ends of the sequence. Must be one character. The default is "!".
    pass_nan : bool, optional
        Whether to pass or raise an exception when given nan in the sequences
    debug : bool, optional
        Whether a message should be printed when failing to find a residue. The default is True.

    Returns
    -------
    list
        [0]: converted residue number
        [1]: alignment score

    """
    
    #print("resnum:", resnum, "\tseq1:", seq1, "\tseq2:", seq2, "\tplaceholder:", placeholder)
    #print(get_residues(seq1[resnum - 1], 5, seq1, placeholder))
    #print(seq1, resnum)
    
    #catch input errors
    assert isinstance(resnum, int), "Expected int for resnum, got " + repr(type(resnum))
    assert isinstance(pass_nan, bool), "Expected bool for pass_nan, got " + repr(type(pass_nan))
    assert isinstance(flanknum, int), "Expected int for flanknum, got " + repr(type(flanknum))
    assert isinstance(placeholder, str), "Expected str for placeholder, got " + repr(type(placeholder))
    assert len(placeholder) == 1, f"Expected single letter string for placeholder, got '{placeholder}'."
    assert isinstance(debug, bool), "Expected bool for debug, got " + repr(type(debug))
    
    if pass_nan == False:
        assert isinstance(seq1, str), "Expected str for seq1, got '" + repr(seq1) + "' which is " + repr(type(seq1))
        assert isinstance(seq2, str), "Expected str for seq2, got '" + repr(seq2) + "' which is " + repr(type(seq2))
    else:
        if isinstance(seq1, str) == False:
            return [np.nan, np.nan]
        if isinstance(seq2, str) == False:
            return [np.nan, np.nan]

    #   Extract the flanking sequence in seq1
    failed = False
    try:
        extract = get_residues(seq1[resnum - 1],
                               flanknum,
                               seq1,
                               placeholder)[resnum]
        if debug == True:
            print("Extracted flanking sequences:", extract)
    except:
        failed = True
        with open("converting regions errors.csv", "a+") as f:
            f.write(f"No residue {resnum} in {seq1} to convert to {seq2}\r")
    
    if failed == False:
        #   Turn seq2 into a dictionary of its relevent residues
        seq2 = get_residues(seq1[resnum - 1],
                            flanknum,
                            seq2,
                            placeholder,
                            frameshift = 1)

        highscore = 0
        highscorer = np.nan
        for potential in seq2:  #   For potential residues in the seq2 dictionary
            score = 0
            for seq2index, letter in enumerate(seq2[potential]):
                group = False
                for key in amino_acid_groups:
                    if letter in amino_acid_groups[key]:
                        group = key
                if letter == extract[seq2index]:
                    score = score + 1
                elif extract[seq2index] in amino_acid_groups[group]:
                    score = score + 0.5
            if debug:
                print("Potential:", potential, "Score:", score)
            if score > highscore:
                highscorer = potential
                highscore = score

        #   highscorer will be a string if it was a frameshifted residue, so turn this
        #   back to a plain integer. Also needs to remove the +/- frameshift.
        if isinstance(highscorer, str):
            try:
                highscorer = int(highscorer.split("-")[0].split("+")[0])
            except ValueError as e:
                print("Failed to convert highscorer to int:", highscorer)
                raise e
        return [highscorer, highscore]
    else:
        return [np.nan, np.nan]
    #except:
    #    if debug == True: print("Failed to get equivalent residue. Resnum:", resnum, " Seq1:", seq1, " Seq2:", seq2)
    #    return [np.nan, np.nan]

def reverse_sequence(sequence, seq_start = False, seq_end = False):
    """
    Takes input characters and reverses it. If seq_start and seq_end are given, function
    returns where seq_start and seq_end are in the new sequence.

    Parameters
    ----------

    sequence : str
        The sequence to reverse.
    seq_start : int, optional
        Region position that will be returned in the new sequence. The default is False.
    seq_end : int, optional
        Region position that will be returned in the new sequence. The default is False.

    Returns
    -------
    str or
    tuple
        [0]: The reversed sequence.
        [1]: The new start position of the region.
        [2]: The new end position of the region.
    """

    # If no start or end is given, just reverse the sequence
    if seq_start is False and seq_end is False:
        return sequence[::-1]
    assert isinstance(seq_start, int), "Expected int for seq_start, got " + repr(type(seq_start))
    assert isinstance(seq_end, int), "Expected int for seq_end, got " + repr(type(seq_end))

    #   If start and end are given, return the new start and end positions
    end = len(sequence)
    new_start = end - seq_end
    new_end = end - seq_start
    return sequence[::-1], new_start + 1, new_end + 1

def convert_region(start_sequence: str,
                   start_region : Union[int, list],
                   end_sequence : str,
                   debug = False) -> dict:
    """
    Takes a start_region in start_sequence and returns where this region is in end_sequence. 
    Use biological sequence numbers (start at 1)

    Parameters
    ----------
    start_sequence : str
        The sequence which the known region belongs to.
    start_region : int | list
        The residue or residue range (as a list) that this region occupies.
    end_sequence : str
        The new sequence where this region should be detected.
    debug : bool, optional
        Should progress be printed.
    
    Returns
    -------
    dict
        Contains "start" and "end" as the start and end of the sequence, and "score" 
        as the alignment score of the new region.

    """

    #   Check if start_region is nan, if so, just return NaN
    if isinstance(start_region, list) is False:
        if pd.isnull(start_region):
            return {"start" : np.nan, "end": np.nan, "score" : np.nan}

    #   Assert some stuff about the sequence
    assert isinstance(start_sequence,
                      str), "Expected string for start_sequence, got " + repr(type(start_sequence))
    if len(start_sequence) == 0:
        raise Exception("start_sequence is empty string.")
    assert type(start_region) in [list, int], "Expected list or int for start_region, got " + repr(type(start_region))
    assert isinstance(end_sequence,
                      str), "Expected string for end_sequence, got " + repr(type(end_sequence))
    if len(end_sequence) == 0:
        raise Exception("end_sequence is emptry string.")

    #   If only one number is presented for start region, duplicate it to give a range
    #   of just one amino acid.
    try:
        #   If list of one, duplicate it.
        if len(start_region) == 1:
            start_region.append(start_region[0])

        #   If one number in the range is np.nan, duplicate the other number into that spot.
        if pd.isnull(start_region[0]):
            start_region[0] = start_region[0]
        if pd.isnull(start_region[1]):
            start_region[1] = start_region[0]
        #   If either number is now np.nan, this means both were np.nan, so just return
        #   nan values.
        if np.nan in start_region:
            return {"start" : np.nan, "end": np.nan, "score" : np.nan}
    except:
        start_region = [int(start_region), (start_region)] # If not a list, duplicate into a list.
    start_sequence = "X" + start_sequence #add an X to the begining of start sequence to fix indexing
    end_sequence = "X" + end_sequence #add an X to the begining of end sequence to fix indexing

    #  Start saving some residues
    residues = []
    try:
        seq_range = range(start_region[0], start_region[1] + 1) #correct for addition of X and make it inclusive
    except:
        raise Exception("""start_region must be number or list of two positive numbers to 
                        specify a sequence region. Got:""" + repr(start_region) + """ for 
                        sequence """ + repr(start_sequence))

    for i in seq_range:
        residue = get_equivalentresidue(resnum = i + 1,
                                        seq1 = start_sequence,
                                        seq2 = end_sequence,
                                        debug = debug)
        residue[0] = residue[0]-1
        if residue[1] > 4:
            residue.append(start_sequence[i])
            residues.append(residue)

    #do this again in reverse
    rev_residues = []
    seq_reversed, new_region_start, new_region_end = reverse_sequence(start_sequence,
                                                                      start_region[0] + 2,
                                                                      start_region[1] + 2)
    reverse_end_sequence = end_sequence[::-1]
    for i in range(new_region_start, new_region_end + 1):
        rev_residue = get_equivalentresidue(resnum = i + 1,
                                            seq1 = seq_reversed,
                                            seq2 = reverse_end_sequence,
                                            debug = debug)
        if rev_residue[1] > 4:
            rev_residue.append(seq_reversed[i])
            rev_residues.append(rev_residue)

    #assign scores for every residue in end_sequence
    end_sequence_scores = [0 for i in end_sequence]
    end_sequence_flanking = [0 for i in end_sequence]

    #and in reverse
    end_sequence_scores_reverse = [np.nan for i in reverse_end_sequence]
    #reverse the reversed scores
    for i in rev_residues:
        i[0] = len(end_sequence_scores) - i[0] #restore the original residue number to each residue
        end_sequence_scores_reverse[i[0]] = i[1] #add scores to every residue in the output
    rev_residues = rev_residues[::-1] #reverses the rev_residues so its in the same order as residues
    #print(residues)

    #combine forward and reverse residues
    for residue in residues:
        end_sequence_scores[residue[0]] = np.nan_to_num(residue[1])
    for residue in rev_residues:
        end_sequence_scores[residue[0]] += np.nan_to_num(residue[1])

    #find the range that the sequence confidently falls under
    try:
        max_score = max(end_sequence_scores)
        first_max_index = False
        last_max_index = 0
        last_max_score = 0
        for index, score in enumerate(end_sequence_scores):
            if score > max_score * 0.7 and first_max_index is False: #begin the definite range where the score is 80% of the max score
                first_max_index = index
                first_max_score = score
            if score > max_score * 0.7: #end the definite range where the score is 80% of the max score
                last_max_index = index
                last_max_score = score
    except:
        raise Exception("sequence range not found.")

    if debug == True:
        print("well defined bounds:", first_max_index, last_max_index)

    # Give residues a score based on all the flanking scores in end_sequence
    for residue in residues:
        flanking_score = 0
        for i in range(-5,6):
            try:
                new_score = end_sequence_scores[residue[0] + i]
                if np.isnan(new_score) is False:
                    flanking_score += new_score
            except:
                pass
        end_sequence_flanking[residue[0]] = flanking_score
        residue.append(flanking_score)

    start = first_max_index
    end = last_max_index

    score = 0
    numbers = 0
    for i in range(start, end + 1):
        if np.isnan(end_sequence_scores[i]) == False:
            score += end_sequence_scores[i]
        numbers += 1
    score = score/numbers
    score = (score/22)*100
    
    return {"start" : start, "end": end, "score" : score}

def int_or_nan(intended_integer):
    """
    
    Parameters
    ----------
    intended_integer : anything convertible to int
        Input that should be converted to int.

    Returns
    -------
    int or float
        Input integer or float if input couldn't be converted.

    """
    try:
        return int(intended_integer)
    except:
        return np.nan

def convert_regions(regions : dict, seq1 : str, seq2) -> dict:
    """
    Takes a dictionary of significant regions as an input, converts them from sequence 1
    to sequence 2. If a dict of sequences is given for seq2, then will find the region 
    in all those sequences.

    Parameters
    ----------
    regions : dict or str in dict format
        Input set of regions to be converted.
    seq1 : str
        Sequence the regions are from.
    seq2 : str or dict
        Sequence the regions should be converted to.

    Returns
    -------
    dict
        New dictionary of regions that have been converted from seq1 sequence to seq2 
        sequence.

    """

    #   Check that none of the inputs are null
    if pd.isnull(regions) or pd.isnull(seq1) or pd.isnull(seq2):
        return np.nan
    if isinstance(regions, str):
        # Convert the dictionary string to dict
        regions = ast.literal_eval(regions) 
    assert isinstance(seq1, str), "Expected str for seq1, got " + repr(type(seq1))
    
    #   Check if seq2 is a dictionary, in which case interpret it as multiple sequences to 
    #   convert to.
    multiple_target_sequences = False
    if isinstance(seq2, dict):
        multiple_target_sequences = True
    else:
        assert isinstance(seq2, str), "Expected str or dict for seq2, got " + repr(type(seq2))
    
    #   If seq2 is only one sequence, wrap it in a dictionary so it is dealt with the 
    #   same as if multiple sequences.
    if multiple_target_sequences == False:
        seq2 = {"placeholder" : seq2}
    
    #   Iterate through the target sequences.
    new_regions_per_seq2 = {}
    
    for target_key in seq2:
        #   Make a dictionary for storing the edited regions
        new_regions = {}
        
        #   Iterate through the regions.
        for region in regions:
            #   Get all the dictionary keys of the regions
            #print("Here we go")
            #print(regions[region])
            keys = regions[region].copy().keys()
            
            #   Copy regions for editing
            new_regions[region] = regions[region].copy()
            
            #   If just a position is offered, indicating 1 amino acid
            if "position" in keys:
                #   Convert that amino acid position
                converted = convert_region(seq1, 
                                           int_or_nan(regions[region]["position"]), 
                                           seq2[target_key]).copy()
                new_regions[region]["position"] = converted.copy()["start"]
                new_regions[region]["score"] = converted.copy()["score"]
                
            #   If a begin and end is offered, so usually an amino acid range
            if "begin" in keys and "end" in keys:
                converted = convert_region(seq1, 
                                           [int_or_nan(regions[region]["begin"]), int_or_nan(regions[region]["end"])],
                                           seq2[target_key])
                new_regions[region]["begin"] = converted.copy()["start"]
                new_regions[region]["end"] = converted.copy()["end"]
                new_regions[region]["score"] = converted.copy()["score"]
                
        #   Save this set of regions in into the dictionary of regions for each seq2
        new_regions_per_seq2[target_key] = new_regions.copy()
    
    #   If seq2 was a string, take the output out of the dictionary it was wrapped in.
    if multiple_target_sequences == False:
        new_regions_per_seq2 = new_regions_per_seq2["placeholder"]
    
    return new_regions_per_seq2

def get_uniprot_sequence(unicode : str,
                         pass_nan : bool = True,
                         pass_no_output = True,
                         debug : bool = False) -> str:
    """
    When given a uniprot accession code will return the sequence of the protein.

    Parameters
    ----------
    unicode : str
        Uniprot accession code.
    pass_nan : bool, optional
        As long as True, input of nan will be passed, otherwise will raise exception.
    pass_no_output : bool, optional
        If false, output of "" will raise exception.
    debug : bool, optional
        Should this function print as it goes.

    Returns
    -------
    str
        Seqeuence of protein to which uniprot accession was given.

    """
    # Check exceptions with the inputs
    for i in [pass_nan, pass_no_output, debug]:
        assert isinstance(i, bool), "get_uniprot_sequence() options must be bools, found " + repr(type(i))
    if pass_nan == False and pd.isnull(unicode):
        raise Exception("Given nan as unicode, which is not allowed if pass_nan == False.")
    if pd.isnull(unicode) == False:
        assert isinstance(unicode, str), "Expected string or nan for unicode, got '" + repr(unicode) + "', which is " + repr(type(unicode))
        
        # Fetch the sequence of the accession given
        if debug == True: print("getting", "https://rest.uniprot.org/uniprotkb/" + unicode + ".fasta")
        soup = requests.get("https://rest.uniprot.org/uniprotkb/" + unicode + ".fasta").text
        output = "".join(soup.split("\n")[1:-1])
        
        # If pass_no_output = False then raise an exception if the output is ""
        if pass_no_output == False and output == "":
            raise Exception(f"No sequence found for {unicode} and pass_no_output is False.")
        
        # Otherwise a message if output is ""
        if output == "" and debug == True:
            print(f"No sequence found for {unicode}.")
        return output
    else: # If input was nan and pass_nan = True then print this message
        if debug == True: print("Given nan as unicode to function get_uniprot_sequence().")

def plusminusorNaN(input, separator):
    """
    This function splits a number with a combined deviation separated (e.g. by ±) into a list 
    with two seperate numbers, or otherwise returns an extra NaN if there is no ±
    """
    try:
        return input.split(separator)
    except:
        return [input, np.nan]

def separatevariance(input, separator, var_addon):
    """
    Takes a pandas dataframe and separates (e.g. by ±) variance into their own columns with a 
    custom header addon.
    """
    data = input
    tosplit = []
    for col in data:
        for bit in data[col]:
            try:
                if separator in bit:
                    tosplit.append(col)
                    break
            except:
                pass
    for col in tosplit:
        data[col + var_addon] = data[col].apply(lambda x : plusminusorNaN(x, separator)[1])
        data[col] = data[col].apply(lambda x : plusminusorNaN(x, separator)[0])
    return data

def get_oximouse_data(age : str):
    """
    Fetch oximouse data as a dataframe. Age can be "aged", "young", or "detected".
    Aged or young will return the data for the aged or young mice, respectively.
    Detected returns a list of every cysteine that was detected with oximouse.
    Oximouse data maps oxidation of specific cysteines in different body parts in
    aged or young mice.
    Reference:
    Xiao, H., Jedrychowski, M. P., Schweppe, D. K., Huttlin, E. L., Yu, Q., Heppner, D. E., Li, J., Long, J., Mills, E. L., Szpyt, J., He, Z., Du, G., Garrity, R., Reddy, A., Vaites, L. P., Paulo, J. A., Zhang, T., Gray, N. S., Gygi, S. P., & Chouchani, E. T. (2020). A Quantitative Tissue-Specific Landscape of Protein Redox Regulation during Aging. Cell, 180(5), 968-983.e24. https://doi.org/10.1016/j.cell.2020.02.012
    """
    if age == "aged":
        df = pd.read_csv("https://piecole.com/data/oximouse_sensitive_old.txt", sep = "\t")
    elif age == "young":
        df = pd.read_csv("https://piecole.com/data/oximouse_sensitive_young.txt", sep = "\t")
    elif age == "detected":
        df = pd.read_csv("https://piecole.com/data/oximouse_all_cysteines_detected.txt", sep = "\t")
    else:
        raise Exception("age must be 'aged', 'young', or 'detected'. Given: " + repr(age) + ".")
    #  Split data into measurement + variance
    df = separatevariance(df, "±", " SEM")

    print(f"Downloaded {age} oximouse data. Cite:")
    print("Xiao, H., Jedrychowski, M. P., Schweppe, D. K., Huttlin, E. L., Yu, Q., Heppner, D. E., Li, J., Long, J., Mills, E. L., Szpyt, J., He, Z., Du, G., Garrity, R., Reddy, A., Vaites, L. P., Paulo, J. A., Zhang, T., Gray, N. S., Gygi, S. P., & Chouchani, E. T. (2020). A Quantitative Tissue-Specific Landscape of Protein Redox Regulation during Aging. Cell, 180(5), 968-983.e24. https://doi.org/10.1016/j.cell.2020.02.012")

    return df

def get_alphafold_structure(uniprot_code : str,
                            folder : str = "structures",
                            extension = "ent",
                            strict = False,
                            debug = False,
                            silent = False):
    """
    Downloads a structure from AlphaFold for a given uniprot code.

    Parameters
    ----------
    uniprot_code : str
        Uniprot code of the protein to download.
    folder : str, optional
        Folder to save the structure to. The default is "structures".
    extension : str, optional
        Extension of the file. The default is "ent".
    strict : bool, optional
        Whether to raise an exception if the structure is not found. The default is False.
    debug : bool, optional
        Whether to print messages as it goes. The default is False.

    Returns
    -------
    None.
    
    """
    folder = parse_folder(folder)
    # Check whether the structure exists
    if os.path.exists(folder + uniprot_code + "." + extension):
        if debug:
            print("Already have structure for " + uniprot_code + ".")
        return

    # Get the structure
    if not silent:
        print("Downloading structure for " + uniprot_code + " from AlphaFold. Please cite: ")
        print("Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2")
    url = "https://alphafold.ebi.ac.uk/files/AF-" + uniprot_code + "-F1-model_v4.pdb"
    data = requests.get(url, allow_redirects=True)

    if "NoSuchKey" in str(data.content):
        if strict:
            raise Exception("AlphaFold structure for " + uniprot_code + " not found.")
        else:
            if debug:
                print("AlphaFold structure for " + uniprot_code + " not found.")
            return

    # Save the structure
    open(folder + uniprot_code + "." + extension, 'wb').write(data.content)

def PDBsearch(query : str) -> list:
    """
    Takes a query and returns a list of PDB IDs that match that query.
    """
    #   Check input is a string
    assert isinstance(query, str), "Expected str for query, got " + repr(type(query))
    
    # Make the query
    if pdb_search_enabled:
        query = PDBquery(query)
    else:
        raise ImportError("PDBsearch() requires the rcsbsearchapi module to be installed.")
    # Execute the query
    results = list(set(query()))

    return results

def get_PDB_structure(pdb_id : str,
                      folder : str = "structures",
                      extension = "ent",
                      debug = False):
    """
    Downloads a structure from the PDB for a given PDB ID.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the structure to download.
    folder : str, optional
        Folder to save the structure to. The default is "structures".
    extension : str, optional
        Extension of the file. The default is "ent".
    strict : bool, optional
        Whether to raise an exception if the structure is not found. The default is False.
    debug : bool, optional
        Whether to print messages as it goes. The default is False.
    Returns
    -------
    None.
    """
    folder = parse_folder(folder)
    # Check whether the structure exists
    if os.path.exists(folder + pdb_id + "." + extension):
        if debug:
            print("Already have structure for " + pdb_id + ".")
        return
    # Get the structure
    if debug:
        print("Downloading structure for " + pdb_id + " from PDB.")
    url = "https://files.rcsb.org/download/" + pdb_id + ".pdb"
    data = requests.get(url, allow_redirects=True)
    open(folder + pdb_id + "." + extension, 'wb').write(data.content)

import propka.run as pk

def run_propka(input_file,
               structure_folder = "pdb",
               structure_extension = "ent",
               propka_folder = "propka/",
               check = True,
               silence = False):
    """
    Checks if a propka file exists, if not then it attempts to make compute one

    Parameters
    ----------

    input_file : str  
        The name of the structure to compute propka for (e.g. "3ii6").
    structure_folder : str, optional  
        The folder to look for the structure in. The default is "pdb".
    structure_extension : str, optional
        The extension of the structure file. The default is "ent".
    propka_folder : str, optional
        The folder to save the propka file in. The default is "propka/".
    check : bool, optional
        Whether to check if the propka file exists before computing it. The default is True.
    silence : bool, optional
        Whether to print messages as it goes. The default is False.
        
    Returns
    -------
    i : propka.run.single
        The propka object. Also saves it to a file.
        
    """
    worked = False

    # Check if propka is installed
    try:
        pk
    except:
        raise Exception("PROPKA not installed. Please install to use this function.")

    propka_folder = parse_folder(propka_folder)
    propka_path = propka_folder + "pdb" + input_file + ".pka"
    # Check if the propka file exists
    if os.path.exists(propka_path) and check:
        if silence == False:
            print(propka_path + f" already exists at {propka_path}. If you want to recompute it, set check = False.")
        return
    else:
        if silence == False:
            print("Going to comput propka for " + input_file + ".")
        structure_folder = parse_folder(structure_folder)
        paths = glob.glob(structure_folder + "/**/" + input_file + "*.*" + structure_extension + "*", recursive = True)
        if paths == []:    
            print("No structure found with glob '" + structure_folder + "/**/" + input_file + "*.*" + structure_extension + "*" + "'.")
            return
        path = paths[0]

        if silence == False:
            print("Calculating pKas with PROPKA and saving to file. Please cite:")
            print("Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values. Chresten R. Søndergaard, Mats H. M. Olsson, Michał Rostkowski, and Jan H. Jensen. Journal of Chemical Theory and Computation 2011 7 (7), 2284-2295. DOI: 10.1021/ct200133y")

        try:
            if "gz" in path:
                with gzip.open(path, "rt") as unzipped:
                    i = pk.single(path.split(structure_extension)[0] + "pdb", optargs = ["-q"], stream = unzipped) #perform PROPKA on the file
            else:
                with open(path, "rt") as f:
                    i = pk.single(path.split(structure_extension)[0] + "pdb", optargs = ["-q"], stream = f)
            worked = True
        except:
            if silence == False:
                print("PROPKA failed for: ", input_file)
            with open("PROPKA failed for.txt", "a") as file:
                file.write(input_file + "\r")
            
    if worked == True:
        # Move the file to propka folder 
        shutil.move(path.split("\\")[-1].split("ent")[0] + "pka", propka_path.replace("/", "/pdb"))
        return i

# CYSTEINE SPECIFIC, MAKE NON CYSTEINE SPECIFIC
def readpropka(filepath): #reads a propka file, saving the cysteines
    found = False
    try:
        with open(filepath, mode = "r") as thisFile:
            found = True
            pkadatareader = csv.reader(thisFile, delimiter = " ")
            pkadata = []
            for pkaline in pkadatareader:
                while '' in pkaline:
                    pkaline.remove('')
                pkadata.append(pkaline)
    except:
        print(filepath, "not found")
    
    if found == True:
        # Save the cysteine rows as a list
        # Previously didn't work with 4-digit residue numbers or more, so changed it.
        cysteines = [i for i in [a for a in pkadata if len(a) > 5] if i[0][0:3] == "CYS"]

        # Deal with the stupid formatting which means some bits of data join onto each other,
        # by splittting these bits of data in two
        for cysteine in cysteines:
            
            # Separate the residue number from the residue name in cases where residue
            # number was digits or more.
            if len(cysteine[0]) > 3:
                cysteine.insert(1, re.split("(\d+)", cysteine[0])[0])
                cysteine.insert(2, re.split("(\d+)", cysteine[0])[1])
                cysteine.pop(0)            

            if len(cysteine) < 16:
                for i in range(7):
                    cysteine.insert(3, '') # add extra spaces to the bond rows
            if len(cysteine) > 13:
                if(len(cysteine[11]) > 3):
                    cysteine.insert(12, re.split("(\d+)", cysteine[11])[0])
                    cysteine.insert(13, re.split("(\d+)", cysteine[11])[1])
                    cysteine.pop(11)
                if(len(cysteine[15]) > 3):
                    cysteine.insert(16, re.split("(\d+)", cysteine[15])[0])
                    cysteine.insert(17, re.split("(\d+)", cysteine[15])[1])
                    cysteine.pop(15)
                if(len(cysteine[19]) > 3):
                    cysteine.insert(20, re.split("(\d+)", cysteine[19])[0])
                    cysteine.insert(21, re.split("(\d+)", cysteine[19])[1])
                    cysteine.pop(19)
        
        #remove the bond rows (only keep the actual cysteine)
        cysteines = [i[:10] for i in cysteines if i[3] != ""]
        
        #turn cysteines into a pandas dataframe
        save_columns = ["resn", "resi", "chain", "pka", "buried", "nil", "desolvation regular 1", "desolvation regular 2", "effects re 1", "effects re 2"]
        data = pd.DataFrame(columns = save_columns, dtype = object)
        for i in cysteines:
            # Turn list into a dataframe column
            newline = pd.DataFrame([i], columns = save_columns)
            data = pd.concat([data, newline],
                             ignore_index = True)
            
        #add the PDBid as a column
        data["PDBid"] = filepath.split("pdb")[-1].split(".pka")[0]
        
        #drop nil
        data = data.drop(columns = ["nil", "resn"])
        
        return data

# MAKING NON CYSTEINE SPECIFIC
def readpropka_better(filepath = None):
    """
    Reads a propka file, extracting PROPKA info about each residue.
    CAN THIS DEAL WITH NMR? CAN PROPKA NOT DEAL WITH NMR?
    """
    found = False
    pkadata = []
    try:
        with open(filepath, mode = "r") as thisFile:
            pkadatareader = csv.reader(thisFile, delimiter = " ")
            for pkaline in pkadatareader:
                while '' in pkaline:
                    pkaline.remove('')
                pkadata.append(pkaline)
            found = True
    except:
        print(filepath, "not found")
    
    if found == True:
        # Save the residues rows as a list
        # Previously didn't work with 4-digit residue numbers or more,
        # so changed it.
        for i in pkadata:
            print(i)
        residues = [i for i in [a for a in pkadata if len(a) > 5] if i[0][0:3] in threeletter and i[0][0:6] != "PROPKA"]
        for i in residues:
            print(i)
        # Deal with the stupid formatting which means some bits of
        # data join onto each other, by splittting these bits of data in two
        for residue in residues:
            # Separate the residue number from the residue name in cases where residue
            # number was digits or more.
            if len(residue[0]) > 3:
                residue.insert(1, re.split("(\d+)", residue[0])[0])
                residue.insert(2, re.split("(\d+)", residue[0])[1])
                residue.pop(0)            

            if len(residue) < 16:
                for i in range(7):
                    residue.insert(3, '') # add extra spaces to the bond rows
            if len(residue) > 13:
                if(len(residue[11]) > 3):
                    residue.insert(12, re.split("(\d+)", residue[11])[0])
                    residue.insert(13, re.split("(\d+)", residue[11])[1])
                    residue.pop(11)
                if(len(residue[15]) > 3):
                    residue.insert(16, re.split("(\d+)", residue[15])[0])
                    print(residue)
                    residue.insert(17, re.split("(\d+)", residue[15])[1])
                    residue.pop(15)
                if(len(residue[19]) > 3):
                    residue.insert(20, re.split("(\d+)", residue[19])[0])
                    residue.insert(21, re.split("(\d+)", residue[19])[1])
                    residue.pop(19)
        
        # Remove the bond rows (only keep the actual residue)
        residues = [i[:10] for i in residues if i[3] != ""]
        
        # Turn list of residues into a pandas dataframe
        save_columns = ["resn",
                        "resi",
                        "chain",
                        "pka",
                        "buried",
                        "nil",
                        "desolvation regular 1",
                        "desolvation regular 2",
                        "effects re 1",
                        "effects re 2"]
        data = pd.DataFrame(columns = save_columns, dtype = object)
        for i in residues:
            newline = pd.DataFrame([i], columns = save_columns)
            data = pd.concat([data, newline],
                             ignore_index = True)
            
        # Add the PDBid as a column
        data["PDBid"] = filepath.split("pdb")[-1].split(".pka")[0]
        
        # Drop nil
        data = data.drop(columns = ["nil"])
        
        return data
    
def read_propkas(folder = "propka/"):
    """
    Read every propka file in a folder and save it as a pandas dataframe.
    
    Parameters
    ----------
    folder : str, optional
        Folder to read propka files from. The default is "propka/".

    Returns
    -------
    pd.DataFrame
    """

    assert isinstance(folder, str), "Expected str for folder, got " + repr(type(folder))

    files = glob.glob(folder + "*.pka")

    datalist = []
    for file in tqdm(files, desc = "Reading propka files"):
        datalist.append(readpropka(file))
    df = pd.concat(datalist, ignore_index = True)

    return df

def check_structure_for_proximal_atoms(structure,
                                       residue_1,
                                       residue_2,
                                       atom_1 = "CA",
                                       atom_2 = "CA",
                                       max_distance = 10,
                                       HSE = False,
                                       b_factor = None,
                                       strict = True,
                                       quiet = True):
    """
    Open a protein structure (or structure) and search for two
    residues that are within a specified distance of each other. Returning a list of
    dictionaries with two residues and their distance from each other.

    Parameters
    ----------
    structure : str or  Bio.PDB.Structure
        Structure or path to the structure file
    residue_1 : str
        Residue three letter code of the first residue
    residue_2 : str
        Residue three letter code of the second residue
    atom_1 : str, optional
        Name of the atom in the first residue. Default is "CA"
    atom_2 : str, optional
        Name of the atom in the second residue. Default is "CA"
    distance : int, optional
        Distance cutoff. Default is 10.
    HSE : bool, optional
        Whether to calculate HSE. Default is False.
    b_factor : string, optional
        Whether to save b-factors. Options are "one", "two". Default is None.
    strict : bool, optional
        Whether to raise an exception if the atoms are not found in
        the intended residues. Default is True.
    quiet : bool, optional
        Whether to print structure parsing errors.

    Returns
    -------
    List  
        List of dictionaries containing the two residues and their distance from each other.

    Examples
    --------
    >>> check_structure_for_proximal_atoms("structures/A2A5R2.ent", "CYS", "CYS", atom_1 = "SG", atom_2 = "SG", max_distance = 5)

    >>> check_structure_for_proximal_atoms("structures/A2A5R2.ent", "CYS", "CYS", atom_1 = "SG", atom_2 = "SG", max_distance = 5, HSE = True)
    
    >>> check_structure_for_proximal_atoms("structures/A2A5R2.ent", "CYS", "CYS", atom_1 = "SG", atom_2 = "SG", max_distance = 5, b_factor = "one")
    """

    assert isinstance(residue_1, str), "Expected str for residue_1, got " + repr(type(residue_1))
    assert isinstance(residue_2, str), "Expected str for residue_2, got " + repr(type(residue_2))
    assert isinstance(atom_1, str), "Expected str for atom_1, got " + repr(type(atom_1))
    assert isinstance(atom_2, str), "Expected str for atom_2, got " + repr(type(atom_2))
    assert isinstance(max_distance, (int, float)), "Expected int or float for max_distance, got " + repr(type(max_distance))
    assert isinstance(HSE, bool), "Expected bool for HSE, got " + repr(type(HSE))
    assert b_factor in [None, "one", "two"], "Expected None, 'one', or 'two' for b_factor, got " + repr(b_factor)

    if isinstance(structure, str):
        # Make sure structure file exists
        assert os.path.exists(structure), "Structure file not found."
        structures = PDBParser(QUIET = quiet).get_structure("structure", structure)
    else:
        # Make sure structure is a structure
        assert isinstance(structure, Structure.Structure), "Expected Bio.PDB.Structure.Structure for structure, got " + repr(type(structure))
        structures = structure
    output_residues = []
    #print("Structures loaded: ", structures)
    for structure in structures:
        # Compile all the residues, so that intermolecular interactions can be
        # checked
        residues = []
        for chain in structure:

            # Save all the b_factors to be added to every residue
            if b_factor is not None:
                b_factors = []
                for residue in chain:
                    # Calculate mean bfactor of atoms in the residue
                    b_factors.append(stats.mean([atom.get_bfactor() for atom in residue]))

            for residue in chain:
                # Only keep relevent residues
                if residue.get_resname() == residue_1 or residue.get_resname() == residue_2:
                    residue.chain = chain
                    
                    if b_factor is not None:
                        residue.b_factor = b_factors

                    residues.append(residue)
        # Iterate through and measure the distance between all the residues
        for residue_A in residues:
            for residue_B in residues:
                # Check the residues are different
                if residue_A != residue_B:
                    # Assert the atoms are in the residues
                    if strict is True:
                        assert atom_1 in [atom.get_id() for atom in residue_A], f"Atom {atom_1} not found in residue {residue_A.get_resname()}{residue_A.get_id()[1]}"
                        assert atom_2 in [atom.get_id() for atom in residue_B], f"Atom {atom_2} not found in residue {residue_B.get_resname()}{residue_B.get_id()[1]}"
                    
                    if atom_1 in [atom.get_id() for atom in residue_A] \
                        and atom_2 in [atom.get_id() for atom in residue_B]:
                        # Measure the distance between the atoms and document if it
                        # is short enough
                        distance = residue_A[atom_1] - residue_B[atom_2]
                        if distance < max_distance:
                            output_dict = {"residue number A" : residue_A.id[1],
                                            "chain A" : residue_A.chain.id,
                                            "residue number B" : residue_B.id[1],
                                            "chain B" : residue_B.chain.id,
                                            "distance" : distance}
                            
                            if HSE is True:
                                res_HSE = get_res_HSE_structure(structure,
                                                                residue_A.chain.id,
                                                                residue_A.id[1],
                                                                residue_B.chain.id,
                                                                residue_B.id[1])
                                
                                output_dict["HSE A num1"] = res_HSE[0][0]
                                output_dict["HSE A num2"] = res_HSE[0][1]
                                output_dict["HSE B num1"] = res_HSE[1][0]
                                output_dict["HSE B num2"] = res_HSE[1][1]

                            if b_factor is not None:
                                if b_factor == "one":
                                    output_dict["b factors"] = residue_A.b_factor
                                if b_factor == "two":
                                    output_dict["b factors A"] = residue_A.b_factor
                                    output_dict["b factors B"] = residue_B.b_factor              

                            output_residues.append(output_dict)
                        
            # Remove residue_A from residues so it wont get tested again
            residues.remove(residue_A)
    return output_residues

def is_amino_acid(residue):
    """
    Check if a residue is an amino acid.
    """
    return residue.get_resname() in threeletter

def combine_range(input : list):
    """
    Take in a series of numbers and output a list of lists where
    adjacent numbers have been combined.
    """
    output = []
    last_number = None
    
    # Convert to list and sort input
    input = list(input)
    input = sorted(input)

    for i in input:
        if last_number is not None:
            if i == last_number + 1:
                output[-1][-1] = i
            else:
                output.append([i,None])
            last_number = i
        else:
            output.append([i,None])
            last_number = i
    for i in output:
        if i[1] is None:
            i.pop(1)
    return output

def extract_interactions(structure,
                         max_distance: int = 4,
                         strict: bool = True,
                         debug=False) -> pd.DataFrame:
    """
    Iterate through the chains in a structure and extract regions that
    interact with ions, ligands, and other chains in the structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Structure to extract interactions from.
    max_distance : int, optional
        Maximum distance for an interaction to be considered. The default is 4.
    strict : bool, optional
        Whether to raise an exception if no interactions are found. If False then
        if no interactions are found will return an empty dictionary. The default is True.

    Returns
    -------
    dict

    Examples
    --------

    >>> from strucseq import strucseq as sq
    >>> structure = sq.PDBParser().get_structure("struc", "structures/3OCP.ent")
    >>> interactions = sq.extract_interactions(structure)

    """
    
    interactions = []

    for model in structure:
        # Collect all residues and atoms in the model
        atoms = []
        atom_to_residue = {}
        for chain in model:
            for residue in chain:
                residue.chain = str(chain.id)
                residue.state = str(model)
                for atom in residue:
                    atoms.append(atom)
                    atom_to_residue[atom] = residue

        if debug:
            print(f"{len(atoms)} atoms found in structure.")
            print(f"Number of atoms collected: {len(atoms)}")
            start_time = time.time()

        # Build the NeighborSearch tree
        if len(atoms) == 0:
            return {}
        neighbor_search = NeighborSearch(atoms)

        # Find all atom pairs within max_distance
        close_atom_pairs = neighbor_search.search_all(max_distance, level='A')

        processed_residue_pairs = set()
        for atom1, atom2 in close_atom_pairs:
            res1 = atom_to_residue[atom1]
            res2 = atom_to_residue[atom2]

            # Skip interactions within the same residue
            if res1 == res2:
                continue

            # Check neither are HOH
            if "HOH" in [res1.get_resname(), res2.get_resname()]:
                continue

            # Skip if residues are in different models (states)
            #if res1.state != res2.state:
            #    continue

            # Skip invalid or duplicate interactions
            if res1 == res2 or res1.state != res2.state \
                or frozenset((res1, res2)) in processed_residue_pairs:
                continue
            processed_residue_pairs.add(frozenset((res1, res2)))

            # Get ligand/ion interactions
            if is_amino_acid(res1) ^ is_amino_acid(res2):
                interactions.extend([{
                    "Chain": res1.chain,
                    "Residue": res1.id[1],
                    "Distance": atom1 - atom2,
                    "Interactor": f"{res2.get_resname()}"
                }, {
                    "Chain": res2.chain,
                    "Residue": res2.id[1],
                    "Distance": atom1 - atom2,
                    "Interactor": f"{res1.get_resname()}"
                }])
            # Get interchain interactions between amino acids
            elif is_amino_acid(res1) and is_amino_acid(res2):
                if res1.chain != res2.chain:
                    interactions.extend([{
                        "Chain": res1.chain,
                        "Residue": res1.id[1],
                        "Distance": atom1 - atom2,
                        "Interactor": f"chain {res2.chain}"
                    }, {
                        "Chain": res2.chain,
                        "Residue": res2.id[1],
                        "Distance": atom1 - atom2,
                        "Interactor": f"chain {res1.chain}"
                    }])

    # Trim the interactions to the shortest interactor atom distance
    df = pd.DataFrame(interactions)
    if df.empty:
        if strict:
            raise ValueError("No intermolecular interactions found in structure. Use strict=False to return an empty dictionary.")
        else:
            return {}
    df = df.sort_values("Distance")
    df = df.drop_duplicates(ignore_index=True, subset=["Chain", "Residue", "Interactor"])

    # Compile interactions into dictionaries with ranges of residues
    df.sort_values(["Interactor", "Chain", "Distance"])
    interactions_dict = {}
    for chain in df["Chain"].unique():
        interactions_dict[chain] = {}
        for interactor in df[df["Chain"] == chain]["Interactor"].unique():
            residues_list = df[(df["Chain"] == chain) & (df["Interactor"] == interactor)]["Residue"].tolist()
            interactions_dict[chain][interactor] = combine_range(residues_list)

    if debug:
        print(f"Interactions extracted in {time.time() - start_time} seconds.")
        # Append to file, create file if needed
        with open("interactions.txt", "a") as file:
            file.write(f"{len(atom_to_residue)}\t{time.time() - start_time}\n")

    return interactions_dict


def get_structure_sequences(structure):
    """
    Takes a structure, returns the sequences of each chain in a dictionary.
    """
    raise DeprecationWarning("Use extract_chain_sequences_from_structure().")
    sequences = {}
    for model in structure:
        for chain in model:
            sequence = ""
            for residue in chain:
                if residue.resname in threetoone:
                    sequence += threetoone[residue.resname]
                else:
                    sequence += "!"
            sequences[chain.id] = sequence
    return sequences

def get_res_HSE_structure(structure,
                            chain1,
                            resn1,
                            chain2 = None,
                            resn2 = None,
                            alternate = False,
                            only_chains = False,
                            fast = True,
                            debug = False):
    model = structure

    if isinstance(model, Structure.Structure):
        model = structure[0]
        
    # Creating feature that speeds up HSE algorithm by removing
    # distant CA atoms before calculating.
    if fast == True:
        save_atoms = []
        # Iterate through every atom, adding them to save_atoms
        # if they are near the relevent one or two.
        for chain in model:
            for res in chain:
                for atom in res:
                    try:
                        if debug:
                            print("Measuring", atom)
                        if atom.id == "CA" and atom - model[chain1][resn1]["CA"] < 13:
                            save_atoms.append(atom)
                        # If two residues provided, check the other
                        # as well
                        if chain2 != None and resn2 != None:
                            if atom.id == "CA" and atom - model[chain2][resn2]["CA"] < 13:
                                save_atoms.append(atom)
                    except: 
                        # Failed to find a target residue
                        #print(f"Failed to use fast HSE.")
                        # Stop trying to use feature that speeds up HSE
                        fast = False
                        if debug:
                            print("Failed fast")
                        break
        
        # Iterate through every CA atom, and remove it if it isn't
        # in save_atoms .
        if fast == True: # Check that saving proximal CAs worked before removing non-proximal CAs   
            for chain in model:
                for res in chain: 
                    try:
                        if res["CA"] not in save_atoms:
                            chain.detach_child(res.id)
                    except:
                        pass
        else:
            pass
            #print("Reverting to computing regular HSE.")
    
    # If we want to ignore other chains in the structure when
    # computing HSE
    if only_chains == True:
        for a in range(len(model) + 1):
            for i in model:
                if i.id not in [chain1, chain2]:
                    model.detach_child(i.id)

    if alternate == False:
        exp_ca = HSExposureCA(model)
        try:
            res_id1 = model[chain1][resn1].get_id()
        except:
            pass
        
        try:
            info1 = exp_ca[(model[chain1].get_id(), res_id1)]
        except:
            info1 = (np.nan, np.nan, 1)
            if debug:
                print("Failed at exp_ca[(model[chain1].get_id(), res_id1)]")
        if chain2 == None:
            return info1
        else:
            try:
                res_id2 = model[chain2][resn2].get_id()
                info2 = exp_ca[(model[chain2].get_id(), res_id2)]
            except:
                info2 = (np.nan, np.nan, 1)
            return [info1, 
                    info2]
    else:
        hse = HSExposure()

def get_res_HSE_file(PDB_file,
                chain1,
                resn1,
                chain2 = None,
                resn2 = None,
                alternate = False,
                only_chains = False,
                fast = True,
                quiet = True): #get the half sphere exposure of the CA atom of a residue: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    
    with gzip.open(PDB_file, "rt") as unzipped:
        structure = PDBParser(QUIET = quiet).get_structure("struc", unzipped)
        return get_res_HSE_structure(structure,
                                chain1,
                                resn1,
                                chain2,
                                resn2,
                                alternate,
                                only_chains,
                                fast)

def assess_intramolecular_confidence(resnum1 : int, resnum2 : int, confidences : list):
    """
    Take two residues (non-pythonic numbering starting at 1) and the list
    of alphafold confidences for residues in the chain and assess the
    confidence of the interaction between the two residues.

    According to alphafold website, anything above 70 is high.

    Parameters
    ----------
    resnum1 : int
        Residue number of the first residue.
    resnum2 : int
        Residue number of the second residue.
    confidences : list
        List of confidences for the chain.

    Returns
    -------
    float
        Confidence of the interaction between the two residues inclusive.

    """

    resnum1 -= 1
    resnum2 -= 1
    lowest = min(resnum1, resnum2)
    highest = max(resnum1, resnum2)

    return min(confidences[lowest - 1:highest])

class Sequence:
    """
    Class for storing and manipulating sequences.
    """

    def __init__(self, sequence = "", **kwargs):
        self.sequence = sequence
        if "sequence_type" in kwargs:
            assert kwargs["sequence_type"] in ["protein", "DNA", "RNA"], "Expected 'protein', 'DNA', or 'RNA' for sequence_type, got " + repr(kwargs["sequence_type"])
            self.sequence_type = kwargs["sequence_type"]
        if "name" in kwargs:
            self.name = kwargs["name"]
    
    # String behaviour of sequence
    def __str__(self):
        return self.sequence
    # Return behaviour
    def __repr__(self):
        return self.sequence
    # Len behaviour
    def __len__(self):
        return len(self.sequence)
    # Subscriptable
    def __getitem__(self, key):
        return self.sequence[key]
    # Adding together behaviour
    def __add__(self, other):
        # Check that the sequence_types are the same
        if hasattr(self, "sequence_type") and hasattr(other, "sequence_type"):
            assert self.sequence_type == other.sequence_type, "Sequence types do not match. " + repr(self.sequence_type) + " and " + repr(other.sequence_type)
        sequence = Sequence(self.sequence + other.sequence)
        if hasattr(self, "sequence_type"):
            sequence.sequence_type = self.sequence_type
        return sequence
    
    # Compair behaviour
    def __eq__(self, other):
        assert isinstance(other, Sequence), "Expected Sequence for other, got " + repr(type(other))
        # If both have sequence types and they do not match, an error will be thrown
        if hasattr(self, "sequence_type") and hasattr(other, "sequence_type"):
            assert self.sequence_type == other.sequence_type, "Sequence types do not match. " + repr(self.sequence_type) + " and " + repr(other.sequence_type)
        return self.sequence == other.sequence

    def reverse(self):
        """
        Reverse the sequence.
        """
        self.sequence = self.sequence[::-1]

    def search_for(self, search_sequence, region = None):
        """
        Search for a sequence within the sequence.
        """
        if region is None:
            region = [1, len(self.sequence)]

        if isinstance(search_sequence, Sequence):
            search_sequence = str(search_sequence)

        print("search_sequence", search_sequence)
        print("region", region)
        print("end_sequence", self.sequence)

        # THIS FUNCTION ASSUMES PROTEIN?
        return convert_region(start_sequence = search_sequence,
                              start_region = region,
                              end_sequence = self.sequence)
    
    def transcribe(self):
        """
        Transcribe DNA to RNA.
        """
        if hasattr(self, "sequence_type"):
            if self.sequence_type == "DNA":
                self.sequence = self.sequence.replace("T", "U")
                self.sequence_type = "RNA"
            else:
                raise Exception("Sequence is not DNA, so cannot be transcribed.")
        else:
            self.sequence = self.sequence.replace("T", "U")
            self.sequence_type = "RNA"

    def blast(self):
        assert hasattr(self, "sequence_type"), "Sequence type must be defined for BLAST search."
        
        if self.sequence_type == "DNA" or self.sequence_type == "RNA":
            result_handle = NCBIWWW.qblast("blastn", "nt", self.sequence)
        elif self.sequence_type == "protein":
            result_handle = NCBIWWW.qblast("blastp", "nr", self.sequence)
        else:
            raise Exception("Sequence type must be 'DNA', 'RNA', or 'protein'. Got " + repr(self.sequence_type))
        return NCBIXML.read(result_handle)

    # Future functions:
        # reverse_complement
        # get_structure (alphafold and or PDB)
        # could put all the sequence functions in this class

def get_human_uniprot_list(filepath = None):
    print("Downloading human uniprot list.")
    df = pd.read_csv("https://piecole.com/data/uniprot_human.txt", sep = "\t")
    if filepath != None:
        assert isinstance(filepath, str), "Expected str for filepath, got " + repr(type(filepath))
        df.to_csv(filepath, sep = "\t", index = False)
    return df