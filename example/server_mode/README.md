### Server mode


To run the server app you first need to create the config.yaml file.
There is an example in this folder:

```yaml
# config.yaml
dataset_files:
  - "../dataset.sdf"
hostname: "0.0.0.0"
port: 8087
api_name: "/roshambo2/search"
verbosity: 2
```

You run it using the roshambo2_server_app program:
```
roshambo2_server_app.py config.yaml
```



Note: this python program should be installed to your path in the install process. If it is not you can run the script with its full path:

```
python <location of roshambo2 git repo>/roshambo2/scripts/roshambo2_server_app.py
```





Once it has started up you can then run searches on the server instances.
Make sure that you set the api url correctly in your search script.
Look at `search_server.py` for an example

```
python search_server.py
```

 In the   `submit_search` function you can provide a dictionary to the options variable.

e.g. 
```
submit_search(url, query_data, options={'optim_mode':'combination', 'max_results':10000, `reduce_over_conformers`:False}, get_structures=True)
```
the parameters will be passed to the normal Roshambo2 compute function



#### Multiple servers

You can run a server app on different nodes. 
Just launch the server program on each machine independently with an appropriate config file.
You can test this on one machine.
In three different terminals run:
```
roshambo2_server_app.py config_node_1.yaml
```
```
roshambo2_server_app.py config_node_2.yaml
```
```
roshambo2_server_app.py config_node_3.yaml
```
You will have 3 servers running with different urls.
You can then search using the `submit_search_multi_server` function.
There is an example in `search_multiple_servers.py`

```
python search_multiple_servers.py
```
The function will send the query to each server and merge the returned results.
