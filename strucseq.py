import time
import requests
from bs4 import BeautifulSoup
import pandas as pd
import lxml
import numpy as np
import gzip
from Bio.PDB import *

try:
    from tqdm import tqdm
    tqdm.pandas()
except:
    def tqdm(iterator, *args, **kwargs):
        return iterator
    
alphabet = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z", "DCY"]
threeletter = ["ALA","DCY","CYS","ASP","GLU","PHE","GLY","HIS","ILE","JXX","LYS","LEU","MET","ASN","OXX","PRO","GLN","ARG","SER","THR","MSE","VAL","TRP","TPO","TYR","SEP"]
threetoone = dict(zip(threeletter, alphabet))

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
    "X":0, #phospohthreonine didn't know what value to give
    "Y":63,
    "Z":0, #phosphothreonine didn't know what value to give
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
def iterate_uniprot_accessions(in_csv : str, chain_cols : list, out_csv : str, delimiter = "\t", rows : str = 10000000000000, debug = True):
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
    
    uniprotpd = uniprotpd.append(previous_file)
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

    
    uniprotpd = uniprotpd.append(previous_file)
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
        print("downloading", "https://rest.uniprot.org/uniprotkb/" + unicode + ".xml")
    except:
        raise Exception("Give a uniprot accession code, recieved " + repr(unicode))

    output = {}
    if debug == True:
        print(f"Extracting Uniprot information from code: {unicode}.")
  
    #   Fetch the xml version of the uniprot site in beautifulsoup
    soup = BeautifulSoup(requests.get("https://rest.uniprot.org/uniprotkb/" + unicode + ".xml").text, "lxml")
    
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
def iterate_uniprot_details(in_csv : str, out_csv : str, uniprot_csv : str = None, species : str = None, debug = True):
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
        
        if len(row["uniprot"]) > 4:
            newinfo = get_uniprot_details(row["uniprot"]) #gets uniprot information as a dictionary
            newinfo["uniprot"] = row["uniprot"] #adds the uniprot accession code to that dictionary
            #creates a 1 row dataframe of that dictionary, with the keys as the columns 
            newinfo = pd.DataFrame([newinfo],  columns = newinfo.keys())
            dataset = dataset.append(newinfo, ignore_index = True) # adds 1 row to the dataset dataframe
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
    
#CHANGE SCREEN TO REFLECT CHANGE
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
            dataset = dataset.append(newinfo, ignore_index = True) # adds 1 row to the dataset dataframe
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
def get_flanking_info(PDB_file : str, debug : bool = False) -> tuple:
    """
        
    Takes a PDB file and returns flanking information for all the cysteines

    Parameters
    ----------
    PDB_file : str
        The path to the input PDB file.
    debug : bool, optional
        Should this function print its output each time.

    Returns
    -------
    tuple
        [0]: Information about flanking residues including sequence, pKas, hydrophobicities, 
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
                    if debug == True: print("res list:",res_list, "residue:", residue)
                    res_list[residue.id[1] -1] = "!" #otherwise add exclamation marks
            for residue in realreslist:
                if residue.get_resname() == "CYS":
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

