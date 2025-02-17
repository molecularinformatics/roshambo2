# Example script to search a roshambo2 server mode instance.
# This example assumes you already have a sever mode program running at:
# "http://localhost:8087/roshambo2/search"
# please look at the user guide for how to do this
# 
# In the main function of this script we specify the URL of the server.
# We then load in a query molecule, convert it to Roshambo2Data format and send it
# via json to the server.


from roshambo2.client import submit_search, prepare_query
from rdkit import Chem

if __name__ == "__main__":

    # Define the URL for the API endpoint
    url = "http://localhost:8086/roshambo2/search"

    
    # load the query from SDF
    query_sdf = "../query.sdf"

    # prepare and convert to roshambo2 format (you should be able to pass an rdkit molecule here instead)
    query_data = prepare_query(query_sdf, color=True)

    # the standard arguments to roshambo2 compute are passed in server mode as an options={} dictionary as key value pairs. For example to set the max
    # number of hits as 10000 and tanimoto_threshold to 0.5.
    options={'optim_mode':'combination', 
        'max_results':10000, 
        'tanimoto_threshold': 0.5}

    results_dict, molecules = submit_search(url, query_data, options=options, get_structures=True)
    
    print(results_dict)

    # molecules is a dict of lists of RDKit molecules. The keys are the name of the query molecules
    # The configurations are the best fit structures, the same as those returned by roshambo2.get_best_fit_structures()

    # We can write them to an SDF file
    qname = list(molecules.keys())[0] # query name is the key of the dict
    mol_list = molecules[qname] # the dict items are lists of molecules

    writer = Chem.SDWriter('best_fit_mols_from_server.sdf')
    for mol in mol_list:
        writer.write(mol)
