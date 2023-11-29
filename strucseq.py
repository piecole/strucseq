import time
import requests
from bs4 import BeautifulSoup
import pandas as pd
import lxml
import numpy as np
import gzip
from Bio.PDB import *
from typing import Union
import ast

try:
    from rcsbsearchapi.search import TextQuery as PDBquery
except:
    print("rcsbsearchapi not installed, some functions may not work.")

try:
    from tqdm import tqdm
    tqdm.pandas()
except:
    def tqdm(iterator, *args, **kwargs):
        return iterator
    print("tqdm not installed, progress bars will not work.")
    
alphabet = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z", "DCY"]
threeletter = ["ALA","DCY","CYS","ASP","GLU","PHE","GLY","HIS","ILE","JXX","LYS","LEU","MET","ASN","OXX","PRO","GLN","ARG","SER","THR","MSE","VAL","TRP","TPO","TYR","SEP"]
threetoone = dict(zip(threeletter, alphabet))

def parse_folder(input_folder : str):
    assert isinstance(input_folder, str), "str expected for input_folder, got" + repr(type(input_folder))
    input_folder = input_folder.replace("\\", "/")
    if input_folder[-1] != "/":
        input_folder = input_folder + "/"
    return input_folder

get_pKa = { #this dictioary defines the pKas of amino acid side chains
    #according to: https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html#prop
    "R": 12.48, "D": 3.65, "C": 8.18, "E": 4.25, "H": 6.00, "K": 10.53, "Y": 10.07, "U":1, "X":5.9, "Z":5.6, "DCY":8.18
    #selenomethionine "U" has been given as pKa 1 as no value is known
    #phosphoserine and phospohthreonine are from 16018962
    #R-cysteine given same value as L-cysteine
}
get_charge = {
    "R": 1, "D": -1, "E": -1, "H": 1, "K": 1
    }

get_hydrophobic7 = { #all the amino acid hydrophobicities at pH7
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
    "U":74, #selenomethonine same as methionine since no value is known for selenomethionine? - check this
    "V":76,
    "W":97,
    "X":0, #phospohthreonine unknown what value to give
    "Y":63,
    "Z":0, #phosphothreonine unknown what value to give
    "DCY":49 #R-cysteine same as L cysteine
}

def get_uniprot_accessions(pdbcode : str, strict = True, selenium = False, debug = False) -> dict:
    """
    Takes a pdb assession code and accesses the PDB to get the Uniprot assession 
    code of the protein for each chain

    Parameters
    ----------
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
    -------
    dict
        Uniprot accession codes of each chain in the strcture. Key is the chain identifier,
        content is the accession code.

    """
    
    #   Catch incorrect inputs as long as strict is True
    if strict == True:
        assert isinstance(pdbcode, str) == True, "Expected 4 letter string for pdbcode, got type: " + type(pdbcode).__name__ + "."
        assert len(pdbcode) == 4, f"Expected 4 letter string for pdbcode, got: '{pdbcode}'."
    
    #   Fetch the xml version of the uniprot site in beautifulsoup
    failed = False
    if selenium == False:
        fails = 0
        worked = False
        while worked == False:
            try:
                if debug == True:
                    print(f"Getting PDB XML for {pdbcode}")
                soup = BeautifulSoup(requests.get("https://files.rcsb.org/view/" + pdbcode + "-noatom.xml").text, "lxml")
                worked = True
            except:
                fails = fails + 1
                time.sleep(fails**2)
        
        entries = {} #  Create a dictionary to store each entry
        failed = False
        try:
            for entry in soup.find("pdbx:struct_ref_seqcategory").find_all("pdbx:struct_ref_seq"): #find each chain entry
                try:
                    info = entry.find("pdbx:pdbx_db_accession").text
                    if len(info) > 4:
                        entries[entry.find("pdbx:pdbx_strand_id").text] = info #extract each one as a uniprot code assigend with a chain letter
                except:
                    if debug == True:
                        print("no uniprot accession found")
        except:
            if debug == True:
                print("no accession section found for ", pdbcode)
            failed = True
            
            for entry in entries:
                if len(entries[entry]) > 6: #   Check if entry is wrong length, then what?
                    failed = True
                    if debug == True:
                        print("Bad uniprots for " + pdbcode)
                    break
            
    return entries


#   OLD VERSION OF THE FUNCTION:
def iterate_uniprot_accessions_OLD(in_csv : str, chain_cols : list, out_csv : str, delimiter = "\t", rows : str = 10000000000000, debug = True):
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
        if debug == True:
            print("already got:", fetched_chains)
    except:
        previous_file = pd.DataFrame()
        fetched_chains = []
        if debug == True:
            print("starting new PDB list")
    
    apd = data[["PDBid", "a chain"]].rename(columns = {"a chain" : "chain"})
    bpd = data[["PDBid", "b chain"]].rename(columns = {"b chain" : "chain"})
    uniquechains = pd.concat([apd, bpd], ignore_index = True).reset_index(drop = True).drop_duplicates()

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
            while tries < 3 and chaindetails == {}: #if no details are returned, try a couple more times
                chaindetails = get_uniprot_accessions(PDBcode) #this is because sometimes no details are returned when there should be some, probably due to anti-scraping measures
                tries = tries + 1
            
            if debug == True:
                if tries == 3:
                    print("no accessions found")
                if 1 < tries < 3:
                    print("found accessions in", str(tries), "tries.")
            if debug == True:
                print(chaindetails)
            
            for letter in chaindetails: #   Iterate each letter in the dictionary
                #   Create a new row with the ID, chain letter, and uniprot code
                uniprotpd.loc[e] = [PDBcode, letter, chaindetails[letter]] 
                e = e + 1
        i = i + 1
    
    uniprotpd = pd.concat([uniprotpd, previous_file])
    try:
        uniprotpd.to_csv(out_csv, sep=delimiter, index = False)
    except:
        if debug == True:
            print("Failed save, using utf-8 instead.")
        uniprotpd.to_csv(out_csv, sep=delimiter, encoding='utf-8', index = False)

def iterate_uniprot_accessions(in_csv : str, chain_cols, out_csv : str, delimiter = "\t", debug = True):
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
        if debug == True:
            print("already got:", fetched_chains)
    except:
        previous_file = pd.DataFrame()
        fetched_chains = []
        if debug == True:
            print("starting new PDB list")
            
    #   Make separate dataframes for each of the chain cols and then combine them
    chain_dfs = []
    uniquechains = pd.DataFrame()
    
    #   Parse chain_cols into a list if it is a string
    if isinstance(chain_cols, str):
        chain_cols = [chain_cols]
    
    for chain_col in chain_cols:
        #   Make the df
        new_df = data[["PDBid", chain_col]].rename(columns = {chain_col : "chain"})
        #   Combine it with the previously made df
        uniquechains = pd.concat([uniquechains, new_df], ignore_index = True).reset_index(drop = True).drop_duplicates()

    uniquePDBs = uniquechains["PDBid"].drop_duplicates().reset_index(drop = True)
    
    uniprotpd = pd.DataFrame({"PDB" : [], "chain": [], "uniprot": []})
    i = 0
    e = 0
    
    pbar = tqdm(np.setdiff1d(list(uniquePDBs), fetched_chains))
    for PDBcode in pbar:
        pbar.set_description("Getting uniprot accessions: " + PDBcode)

        chaindetails = {}
        tries = 0
        while tries < 3 and chaindetails == {}: #if no details are returned, try a couple more times
            chaindetails = get_uniprot_accessions(PDBcode) #this is because sometimes no details are returned when there should be some, probably due to anti-scraping measures
            tries = tries + 1
        
        if debug == True:
            if tries == 3:
                print("no accessions found")
            if 1 < tries < 3:
                print("found accessions in", str(tries), "tries.")
        if debug == True:
            print(chaindetails)
        
        for letter in chaindetails: #   Iterate each letter in the dictionary
            #   Create a new row with the ID, chain letter, and uniprot code
            uniprotpd.loc[e] = [PDBcode, letter, chaindetails[letter]] 
            e = e + 1

    
    uniprotpd = pd.concat([uniprotpd, previous_file])
    try:
        uniprotpd.to_csv(out_csv, sep=delimiter, index = False)
    except:
        if debug == True:
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
        if debug == True:
            print("downloading", url)
    except:
        raise Exception("Give a uniprot accession code, recieved " + repr(unicode))

    output = {}
    if debug == True:
        print(f"Extracting Uniprot information from code: {unicode}.")
  
    #   Fetch the xml version of the uniprot site in beautifulsoup
    soup = BeautifulSoup(requests.get(url).text, "lxml")
    
    #   Get the protein name.
    try:
        output["uniprot name"] = soup.find_all("fullname")[0].text
    except:
        output["uniprot name"] = np.NaN
    
    #   Get the uniprot abbreviation.
    try:
        output["uniprot abbreviation"] = soup.find_all("name")[0].text
    except:
        output["uniprot abbreviation"] = np.NaN
    
    try:
        #   Search for the list of localisations and then put each one in a list
        localisations = [x.text for x in soup.find("comment", {"type" : "subcellular location"}).find_all("location")]
        #   Output this list
        output["localisations"] = localisations
    except:
        #   If there is no such localisation
        output["localisations"] = np.NaN
    
    #   Get the first description
    try:
        output["description"] = soup.find("title").text
    except:
        #   If no description found.
        output["description"] = np.NaN
        
    #   Get disease variants
    variants = {}
    for variant in soup.find_all("feature", {"type":"sequence variant"}):
        try:
            variant["description"] # Check the variant has a description
        except:
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
                    test = region.find(location)["position"]
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
     
    if debug == True:
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
        if uniprot_csv != None: #   If uniprotdata was specified, then add this to the dataframe
            data = pd.merge(data, uniprotdata, how = "left", left_on=["PDBid", chain + " chain"], right_on=["PDB", "chain"]) #  Looks at the chain letter and pdb accession code and adds the row from the uniprot dataframe
            data = data.drop(["PDB", "chain"], axis = 1) #  Removes the added PDB code and chain since this information is duplicate
            if species == None:
                #   Renames the new 'chain' column to specify a or b
                data = data.rename(columns = {"uniprot" : "uniprot " + chain}) #    
            else:
                data = data.rename(columns = {"uniprot" : chain + " " + species + " accession code"})
    
    data = data.loc[:,~data.columns.duplicated()].copy()
    
    #if a species is specified, then look in that column for the accession code otherwise look in default column
    accession_columns = ["a", "b"]
    if species == None:
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
            if debug == True:
                print("not a uniprot accession:", row["uniprot"])
            
            
    #convert species column into dictionary for better access
    accession_columns = {"a chain": accession_columns[0], "b chain" : accession_columns[1]}
    if species == None: #convert the species into text that can be used to name columns
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
    
def iterate_uniprot_details(in_csv : str, chain_cols : list, out_csv : str, delimiter : str = "\t", uniprot_csv : str = None, species : str = None, debug = True):
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
    for chain_col in chains_cols: #  Adds the uniprot based on on 'a chain' then 'b chain'
        if uniprot_csv != None: #   If uniprot_csv was specified, then add this to the dataframe
            data = pd.merge(data, uniprot_csv, how = "left", left_on=["PDBid", chain_col], right_on=["PDB", "chain"]) #  Looks at the chain letter and pdb accession code and adds the row from the uniprot dataframe
            data = data.drop(["PDB", "chain"], axis = 1) #  Removes the added PDB code and chain since this information is duplicate
            if species == None:
                #Renames the new 'chain' column to specify the chain
                data = data.rename(columns = {"uniprot" : "uniprot " + chain_col}) #    
            else:
                data = data.rename(columns = {"uniprot" : chain_col + " " + species + " accession code"})
    
    data = data.loc[:,~data.columns.duplicated()].copy()
    
    #if a species is specified, then look in that column for the accession code otherwise look in default column
    accession_columns = ["a", "b"]
    if species == None:
        accession_columns = ["uniprot " + x for x in chains] #default column
    else:
        try:
            accession_columns = [x + " " + species + " accession code" for x in chains]
        except:
            raise Exception("Expected species as string, got " + repr(species))
    
    #accesses the internet using the new uniprots to get uniprot data using get_uniprot_details()
    #creates a new series of unique uniprot codes
            
    #GOT TO ABOUT HERE FOR EDITING
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
        
        if len(row["uniprot"]) > 4:
            newinfo = get_uniprot_details(row["uniprot"]) #gets uniprot information as a dictionary
            newinfo["uniprot"] = row["uniprot"] #adds the uniprot accession code to that dictionary
            #creates a 1 row dataframe of that dictionary, with the keys as the columns 
            newinfo = pd.DataFrame([newinfo],  columns = newinfo.keys())
            dataset = pd.concat([dataset, newinfo], ignore_index = True) # adds 1 row to the dataset dataframe
        else:
            if debug == True:
                print("not a uniprot accession:", row["uniprot"])
            
            
    #convert species column into dictionary for better access
    accession_columns = {"a chain": accession_columns[0], "b chain" : accession_columns[1]}
    if species == None: #convert the species into text that can be used to name columns
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

parser = PDBParser()
def get_flanking_info(PDB_file : str, amino_acid : str, debug : bool = False) -> tuple:
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
    with gzip.open(PDB_file, "rt") as unzipped: #open the structure
        structure = parser.get_structure("struc", unzipped) #parse the structure
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
                    if debug == True: print("res list:", res_list, "residue:", residue)
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

def get_residues(residue : str, flanknum : int, sequence : str, placeholder : str = "!", frameshift : bool = False) -> dict:
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

    Returns
    -------
    dict
        A dictionary containing sequence numbers as keys with the data attached to each 
        key being a sequence around that residue.
        If frameshift is True, some entries will be frameshifted sequences with a key
        as the residue number followed by +/- and a number denoting the number of 
        frameshifts.

    """

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
        for res in reslist:    #   For every seq1 residue
        
            #   Set the frameshift level to the flanknum if frameshift not specified
            if isinstance(frameshift, bool):
                frameshift = flanknum
            
            #   Frameshift every potential residue on the left side
            #   For each number of frameshifts up to the flanknum
            for shift in range(frameshift):
                #   Make a new entry to store the shifted sequence as a list
                shifted[str(res) + "-" + str(shift)] = list(reslist[res])
                #   Add an X for every frameshift and remove a residue to compensate
                for i in range(shift):
                    shifted[str(res) + "-" + str(shift)].insert(flanknum, "X")
                    shifted[str(res) + "-" + str(shift)].pop(0)
                    
            #   Frameshift every potential residue on the right side
            for shift in range(frameshift):
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
def get_equivalentresidue(resnum : int, seq1 : str, seq2 : str, flanknum : int = 5, placeholder : str = "!", pass_nan : bool = True, debug : bool = False) -> list:
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
    if pass_nan == False:
        assert isinstance(seq1, str), "Expected str for seq1, got '" + repr(seq1) + "' which is " + repr(type(seq1))
        assert isinstance(seq2, str), "Expected str for seq2, got '" + repr(seq2) + "' which is " + repr(type(seq2))
    assert isinstance(flanknum, int), "Expected int for flanknum, got " + repr(type(flanknum))
    assert isinstance(placeholder, str), "Expected str for placeholder, got " + repr(type(placeholder))
    assert len(placeholder) == 1, f"Expected single letter string for placeholder, got '{placeholder}'."
    assert isinstance(debug, bool), "Expected bool for debug, got " + repr(type(debug))
    
    #   Extract the flanking sequence in seq1
    failed = False
    try:
        extract = get_residues(seq1[resnum - 1], flanknum, seq1, placeholder)[resnum]
    except:
        failed = True
        with open("converting regions errors.csv", "a+") as f:
            f.write(f"No residue {resnum} in {seq1} to convert to {seq2}\r")
    
    if failed == False:
        #   Turn seq2 into a dictionary of its relevent residues
        seq2 = get_residues(seq1[resnum - 1], flanknum, seq2, placeholder, frameshift = 1)
        
        highscore = 0
        highscorer = np.NaN
        for potential in seq2:  #   For potential residues in the seq2 dictionary
            score = 0
            for seq2index, letter in enumerate(seq2[potential]):
                group = False
                for key in amino_acid_groups:
                    if letter in amino_acid_groups[key]:
                        group = key
                if(letter == extract[seq2index]):
                    score = score + 1
                elif extract[seq2index] in amino_acid_groups[group]:
                    score = score + 0.5
                    
            if score > highscore:
                highscorer = potential
                highscore = score
                
        #   highscorer will be a string if it was a frameshifted residue, so turn this
        #   back to a plain integer. Also needs to remove the +/- frameshift.
        if isinstance(highscorer, str):
            highscorer = int(highscorer.split("-")[0].split("+")[0])
                
        return [highscorer, highscore]
    else:
        return [np.NaN, np.NaN]
    #except:
    #    if debug == True: print("Failed to get equivalent residue. Resnum:", resnum, " Seq1:", seq1, " Seq2:", seq2)
    #    return [np.NaN, np.NaN]

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
    if seq_start == False and seq_end == False:
        return sequence[::-1]
    assert isinstance(seq_start, int), "Expected int for seq_start, got " + repr(type(seq_start))
    assert isinstance(seq_end, int), "Expected int for seq_end, got " + repr(type(seq_end))

    #   If start and end are given, return the new start and end positions
    end = len(sequence)
    new_start = end - seq_end
    new_end = end - seq_start
    return sequence[::-1], new_start + 1, new_end + 1

def convert_region(start_sequence: str, start_region : Union[int, list], end_sequence : str, debug = False) -> dict:
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
    if isinstance(start_region, list) == False:
        if pd.isnull(start_region):
            return {"start" : np.NaN, "end": np.NaN, "score" : np.NaN}
    
    #   Assert some stuff about the sequence
    assert isinstance(start_sequence, str), "Expected string for start_sequence, got " + repr(type(start_sequence))
    assert type(start_region) in [list, int], "Expected list or int for start_region, got " + repr(type(start_region))
    assert isinstance(end_sequence, str), "Expected string for end_sequence, got " + repr(type(end_sequence))
    
    #   If only one number is presented for start region, duplicate it to give a range 
    #   of just one amino acid.
    try:
        #   If list of one, duplicate it.
        if len(start_region) == 1:  
            start_region.append(start_region[0])      
            
        #   If one number in the range is np.NaN, duplicate the other number into that spot.
        if pd.isnull(start_region[0]):
            start_region[0] = start_region[0]
        if pd.isnull(start_region[1]):
            start_region[1] = start_region[0]
        #   If either number is now np.NaN, this means both were np.NaN, so just return 
        #   nan values.
        if np.NaN in start_region:
            return {"start" : np.NaN, "end": np.NaN, "score" : np.NaN}
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
        residue = get_equivalentresidue(i + 1, start_sequence, end_sequence)
        residue[0] = residue[0]-1
        if residue[1] > 4:
            residue.append(start_sequence[i])
            residues.append(residue)
    
    #do this again in reverse
    rev_residues = []
    seq_reversed, new_region_start, new_region_end = reverse_sequence(start_sequence, start_region[0] + 2, start_region[1] + 2)
    reverse_end_sequence = end_sequence[::-1]
    for i in range(new_region_start, new_region_end + 1):
        rev_residue = get_equivalentresidue(i + 1, seq_reversed, reverse_end_sequence)
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
            if score > max_score * 0.7 and first_max_index == False: #begin the definite range where the score is 80% of the max score
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
                if np.isnan(new_score) == False:
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
        return np.NaN

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
        return np.NaN
    if isinstance(regions, str):
        regions = ast.literal_eval(regions) # Convert the dictionary string to dict
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

def get_uniprot_sequence(unicode : str, pass_nan : bool = True, pass_no_output = True, debug : bool = True) -> str:
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
        return [input, np.NaN]

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
        raise Exception("age must be 'aged', 'young', or 'detected'.")
    #  Split data into measurement + variance
    df = separatevariance(df, "±", " dev")

    print(f"Downloaded {age} oximouse data. Please cite:")
    print("Xiao, H., Jedrychowski, M. P., Schweppe, D. K., Huttlin, E. L., Yu, Q., Heppner, D. E., Li, J., Long, J., Mills, E. L., Szpyt, J., He, Z., Du, G., Garrity, R., Reddy, A., Vaites, L. P., Paulo, J. A., Zhang, T., Gray, N. S., Gygi, S. P., & Chouchani, E. T. (2020). A Quantitative Tissue-Specific Landscape of Protein Redox Regulation during Aging. Cell, 180(5), 968-983.e24. https://doi.org/10.1016/j.cell.2020.02.012")

    return df

def get_alphafold_structure(uniprot_code, folder):
    # Get the structure
    url = "https://alphafold.ebi.ac.uk/files/AF-" + uniprot_code + "-F1-model_v4.pdb"
    data = requests.get(url, allow_redirects=True)
    # Save the structure
    open(parse_folder(folder) + uniprot_code + ".pdb", 'wb').write(data.content)

def PDBsearch(query : str) -> list:
    """
    Takes a query and returns a list of PDB IDs that match that query.
    """
    #   Check input is a string
    assert isinstance(query, str), "Expected str for query, got " + repr(type(query))
    
    # Make the query
    query = PDBquery(query)
    # Execute the query
    results = set(query())

    return results