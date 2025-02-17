import requests
import json
import os
import sys
from roshambo2.classes import Roshambo2DataReaderSDF, Roshambo2DataReaderRDKit, Roshambo2DataReaderh5
import pickle
import base64
import pandas as pd
import time
from rdkit.Chem.rdchem import Mol as rdkitMol
from roshambo2.pharmacophore import PharmacophoreGenerator

# TODO: change print to logging

def prepare_query(queries, color=False, color_generator=None,  remove_Hs_before_color_assignment=False, n_cpus_prepare=1):
    """Prepare query molecule(s) for a Roshambo2Server search.

        Args:
            queries (str | rdkitMol | List[str] | List[rdkitMol]): The query molecules, can be an SDF filename, Roshambo2H5 filename, an RDKit molecule, or a list of one of these types.
            color (bool, Optional): If True color features will be assinged using the color_generator.
            color_generator (PharmacophoreGenerator, Optional). If passed this color_generator will be used. If None (the default) the Roshambo2 default color generator will be used. *Note the color_generator here must have the same features as the one used by the Roshambo2Server.*

        Returns:
            Roshambo2Data: An Roshambo2Data object.
    
    """

    # loads an SDF or an RDKit mol

    # if user wants color but does not provide a generator we use the default one
    if color_generator is None:
        if color == True:
            color_generator = PharmacophoreGenerator()
        else: 
            color_generator = None
    else:
        color_generator = color_generator
    #TODO: server/database color definitions must match these ones
    #TODO: interaction matrix can be different and should be able to be set here and send with the query

    if isinstance(queries, list):
        if all(isinstance(query, rdkitMol) for query in queries):
            query_data_reader = Roshambo2DataReaderRDKit(queries, color_generator=color_generator, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
        elif all( query.endswith('.sdf') for query in queries):
            query_data_reader = Roshambo2DataReaderSDF(queries, color_generator=color_generator, conformers_have_unique_names=True, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
        elif all(query.endswith('.h5') for query in queries):
            query_data_reader = Roshambo2DataReaderh5(queries)

        else:
            raise ValueError('query input list must be all ".sdf" files or all ".h5" files, or all RDKit molecules')
            
    else:
        if isinstance(queries, rdkitMol):
            query_data_reader = Roshambo2DataReaderRDKit([queries], color_generator=color_generator, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
        elif queries.endswith('.sdf'):
            query_data_reader = Roshambo2DataReaderSDF([queries], color_generator=color_generator, conformers_have_unique_names=True, keep_original_coords=True, remove_Hs_before_color_assignment=remove_Hs_before_color_assignment, n_cpus=n_cpus_prepare)
        elif queries.endswith('.h5'):
            query_data_reader = Roshambo2DataReaderh5([queries])
        
        else:
            raise ValueError('query file must be ".sdf", ".h5", or an RDKit molecule')
    
      

    # merge queries into one dataset
    queries = [query for query in query_data_reader.get_data()]
    query = queries[0]
    for i in range(1,len(queries)):
        query+=queries[i]

    # the queries must all have unique names for the output dict format to work correctly.
    seen  = set()
    for qname in query.f_names:
        if qname in seen:
            raise ValueError("query molecules must have unique names")
        else:
            seen.add(qname)
    

    return query


def check_status(server_url, task_id):
    """Check the status of a running task.

        Args:
            sever_url (str): url of the running server.
            task_id (str): ID of the task.

        Returns:
            str: task status, will return None if the status cannot be got.
    """
    status_url = f'{server_url}/status/{task_id}'
    response = requests.get(status_url)

    if response.status_code == 200:
        return response.json()['status']
    else:
        print(f"Error checking status: {response.json()}")
        return None

def get_result(server_url, task_id, get_structures):
    """Get the result of a task from the server.
    
        Args:
            sever_url (str): url of the running server.
            task_id (str): ID of the task.
            get_structures (bool): If True the RDKit molecules of the best fit results will be returned along with the scores.
    
        Returns:
            Tuple(Dict,Dict): Returns two dictionaries. One with the Pandas dataframe of scores for each query and the other with a list of best fit molecules (as RDKit molecules) for each query.

    """

    if get_structures:
        result_url = f'{server_url}/results_with_structures/{task_id}'
        response = requests.get(result_url)

        if response.status_code == 200:
            
            # get data as json
            data = response.json()
            
            # parse the molecules
            molecules_dict = {}
            for key, encoded_binary_mol_list in data['molecules'].items():
                molecules_dict[key] = []

                for eb_mol in encoded_binary_mol_list:

                    # decode
                    b_mol = base64.b64decode(eb_mol)

                    # from rdkit binary to rdkit
                    mol = rdkitMol(b_mol)

                    # add server url property
                    mol.SetProp('server_url', server_url)

                    molecules_dict[key].append(mol)


            # parse the dataframe
            dataframes_dict = {}
            for key, json_str in data['scores'].items():
                
                dataframe = pd.read_json(json_str) # TODO: FutureWarning: Passing literal json to 'read_json' is deprecated and will be removed in a future version. To read from a literal string, wrap it in a 'StringIO' object.

                # add the server where the data came from:
                dataframe['server_url'] = server_url

                dataframes_dict[key] = dataframe

            return dataframes_dict, molecules_dict
        else:
            print(f"Error retrieving result: {response}")
            return None, None


    else: # just get the scores dict:
        result_url = f'{server_url}/result/{task_id}'
        response = requests.get(result_url)

        if response.status_code == 200:
            # convert results from json back to pandas
            dataframes_dict = {}
            for key, json_str in response.json().items():
               
                dataframe = pd.read_json(json_str) # TODO: FutureWarning: Passing literal json to 'read_json' is deprecated and will be removed in a future version. To read from a literal string, wrap it in a 'StringIO' object.

                # add the server where the data came from:
                dataframe['server_url'] = server_url

                dataframes_dict[key] = dataframe


            return dataframes_dict,None
        else:
            print(f"Error retrieving result: {response.json()}")
            return None, None

def get_error(server_url, task_id):
    """Get the error of a failed task.

        Args:
            sever_url (str): url of the running server.
            task_id (str): ID of the task.

        Returns:
            str: Error message, will return None if there is an error getting the error message.
    """

    result_url = f'{server_url}/result/{task_id}'
    response = requests.get(result_url)

    if response.status_code == 200:
        return str(response.json())
    else:
        print(f"Error retrieving result: {response.json()}")
        return None

def submit_search(server_url, query_data, options=None, get_structures=True):
    """Submit a search on an Roshambo2Server instance.

        Args:
            server_url (str): url of the running server.
            query_data (Roshambo2Data): prepared query data.
            options (Dict, Optional): Dictionary of keyword arguments to send to the Roshambo2.compute() method.
            get_structures (bool, Optional): If True (the default) the best fit molecules will be returned as RDkit molecules. If False only the scores are returned.

        Returns:
            Tuple(Dict, Dict): Returns two dictionaries. One with the Pandas dataframe of scores for each query and the other with a list of best fit molecules (as RDKit molecules) for each query.

    """
    #TODO: probably need some timeout commands in all the post or get requests in-case connection drops


    print("Submitting search...")
    if options is None:
        options = {}

    else: #TODO more validation checks
        for option in ['backend','n_gpus', 'write_scores']:
            if option in options:
                print("option {option} is ignored in server mode and is controlled on the server side.")
                del options[option]    

    # Serialize the query_data with pickle and encode it with base64
    pickled_data = pickle.dumps(query_data)
    encoded_query_data = base64.b64encode(pickled_data).decode('utf-8')

    payload = {
        'query_data': encoded_query_data,
        'options': options,
        'get_structures': get_structures
    }

    # Send the POST request to submit the search
    response = requests.post(server_url, json=payload)
    
    time.sleep(1)


    if response.status_code == 202:
        task_id = response.json()['task_id']
    else:
        print(f"Error submitting search: {response.json()}")
        return None
    
    if task_id:
        print(f"Search submitted successfully. Task ID: {task_id}")
        
        # poll status
        while True:
            status = check_status(server_url, task_id)
            if status == 'completed':
                print(f"Task {task_id} completed. Fetching results...")
                sucess=True
                break
            elif status == 'failed':
                print(f"Task {task_id} failed.")
                failed_result = get_error(server_url, task_id)
                print("Error from server: ", failed_result)
                sucess=False
                break
            else:
                print(f"Task {task_id} is still {status}. Checking again in 5 seconds...")
                time.sleep(5)
        
        if sucess:
            result = get_result(server_url, task_id, get_structures=get_structures)
            if result:
                return result
            else:
                print("Failed to retrieve result")
                return None
        else:
            return None



def _merge_server_results(results, max_results):
    """merges results from multi server search into one via reduction"""


    # results is a list of tuples
    # tuple[0] is a dict of str:pandas dataframe, tuple[1] is a dict of str: list of rdkit molecules or None.


    # we need to merge for each query

    # first result list, tuple 0
    qnames = results[0][0].keys()
    print(qnames)

    output_dfs={}
    output_molecules={}

    # loop over queries
    for qname in qnames:
        # loop over results
    
        dfs = []
        molecules_lists = []

        for (df,molecules_list) in results:

            dfs.append(   df[qname]  )
            molecules_lists.append(   molecules_list[qname]  )

            

        # use pandas to merge
        merged_df = pd.concat(dfs)
        merged_molecules = sum(molecules_lists,[])
        
        # sort by score
        sorted_indices = merged_df['tanimoto_combination'].argsort()[::-1]
        sorted_df = merged_df.iloc[sorted_indices]
        sorted_df = sorted_df.reset_index(drop=True)
        sorted_molecules = [merged_molecules[i] for i in sorted_indices]  # Sort the list using the same index order


        # take max
        if max_results is not None:
            sorted_df = sorted_df.head(max_results)
            sorted_molecules = sorted_molecules[:max_results]
        
        # add to per query dicts
        output_dfs[qname] = sorted_df
        output_molecules[qname] = sorted_molecules

    return output_dfs, output_molecules


def submit_search_multi_server(server_urls, query_data, options=None, get_structures=True):
    """Submit a search on an Roshambo2Server instance.

        Args:
            server_urls (list): list of urls of the running servers.
            query_data (Roshambo2Data): prepared query data.
            options (Dict, Optional): Dictionary of keyword arguments to send to the Roshambo2.compute() method.
            get_structures (bool, Optional): If True (the default) the best fit molecules will be returned as RDkit molecules. If False only the scores are returned.

        Returns:
            Tuple(Dict, Dict): Returns two dictionaries. One with the Pandas dataframe of scores for each query and the other with a list of best fit molecules (as RDKit molecules) for each query.

    """
    #TODO: probably need some timeout commands in all the post or get requests incase connection drops


    print("Submitting search...")
    if options is None:
        options = {}

    else: #TODO more validation checks
        for option in ['backend','n_gpus', 'write_scores']:
            if option in options:
                print("option {option} is ignored in server mode and is controlled on the server side.")
                del options[option]    

    # remember maxn
    if "max_results" in options:
        max_results  = options["max_results"]
    else:
        max_results = None

    # Serialize the query_data with pickle and encode it with base64
    pickled_data = pickle.dumps(query_data)
    encoded_query_data = base64.b64encode(pickled_data).decode('utf-8')

    payload = {
        'query_data': encoded_query_data,
        'options': options,
        'get_structures': get_structures
    }


    responses=[]
    for server_url in server_urls:
        # Send the POST request to submit the search
        response = requests.post(server_url, json=payload)
        responses.append(response)
    
    
    time.sleep(1)

    task_ids=[]
    for response in responses:

        if response.status_code == 202:
            task_id = response.json()['task_id']

            if task_id:
                print(f"Search submitted successfully. Task ID: {task_id}")
                task_ids.append(task_id)

            else:
                return None
        else:
            print(f"Error submitting search: {response.json()}")
            return None
    
    
    completed = [False for _ in task_ids]
    results=[None for _ in task_ids]
    
    # poll status
    while not all(completed):
        for i,(server_url,task_id) in enumerate(zip(server_urls, task_ids)):
            status = check_status(server_url, task_id)
            if status == 'completed':
                print(f"Task {task_id} on server {server_url} completed. Fetching results...")
                
                result = get_result(server_url, task_id, get_structures=get_structures)
                if result:
                    results[i]=result
                    completed[i]=True
                else:
                    print(f"Failed to retrieve result from server {server_url}")
                    completed[i]=True


            elif status == 'failed':
                print(f"Task {task_id} on server {server_url} failed.")
                failed_result = get_error(server_url, task_id)
                print("Error from server: ", failed_result)
                completed[i]=True
            
            else:
                print(f"Task {task_id} on server {server_url} is still {status}.")
        

        print(f"Checking again in 5 seconds...")
        time.sleep(5)

    # if we are here all must be completed, but not necessarily successfully
    if None in results:
        print("One of the searches failed")
        return None
    else:

        return _merge_server_results(results, max_results)