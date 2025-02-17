# User Guide

## Introduction
Roshambo2 is a python package for 3D molecular shape comparison. It uses the Gaussian overlap method of Grant and Pickup, 1995, https://doi.org/10.1021/j100011a016.

It is an evolution of the ROSHAMBO software (https://pubs.acs.org/doi/abs/10.1021/acs.jcim.4c01225, https://github.com/molecularinformatics/roshambo). The core shape overlap optimizer has been re-implemented in C++ and CUDA and the python wrapper code has been re-written. The emphasis is on high speed 3D similarity scoring on large scale datasets.

To facilitate the high performance, and ease of integration with existing chemoinformatics workflows, the program is implemented as a Python Package. You will need to write Python scripts to run your own scoring. We provide simple examples to get you started.


## Getting started
First you will need to install the software. Currently this must be done from source, please see the Install section of the documentation.

In future we will provide PyPI and conda-forge binary distributions.

## Basic usage

In a typical use case, you provide one or more query molecules along with a larger set of dataset molecules. The software then computes the 3D similarity between them, ranking the dataset molecules based on their similarity.

Roshambo2 does not perform 3D conformer generation. You will need to generate the conformers yourself beforehand. 

The basic input to Roshambo2 are 3D SDF files. One SDF for the query molecule(s) and another for the dataset molecules. 

To run the scoring you then need to make a simple python script and run it.

```python
from roshambo2 import Roshambo2
 
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf")
scores = roshambo2_calculator.compute() 
```
The above script is in the example folder of this repo
```
>>> git clone <roshambo2.git>
>>> cd roshambo2/example 
>>> python simple_run.py
```
You will see output like this:

```
INFO:roshambo2.prepare:input file = dataset.sdf, out file = None, assigning color features = False
INFO:roshambo2.prepare:Reading dataset.sdf which contains 153 molecules
INFO:roshambo2.prepare:input file = query.sdf, out file = None, assigning color features = False
INFO:roshambo2.prepare:Reading query.sdf which contains 1 molecules
INFO:roshambo2.roshambo2:Roshambo2 setup completed in 0.26940667105372995s
INFO:roshambo2.roshambo2:Starting compute
INFO:roshambo2.roshambo2:Initializing backend
INFO:roshambo2.roshambo2:Staring optimization
INFO:roshambo2.roshambo2:completed optim of 1 x 153 molecules in 0.12250700802542269 seconds
INFO:roshambo2.roshambo2:saving scores for query CHEMBL221029 to scores_query_CHEMBL221029.csv
INFO:roshambo2.roshambo2:Roshambo2 compute completed in 0.141 seconds. Searched 153 configs.
```

It will have saved the scores to the csv file named `scores_query_CHEMBL221029.csv`. The name `CHEMBL221029` in this example comes from the name field in the `query.sdf` file. If there were more molecules in the query file you would get an output scores csv for each one named correspondingly.

### More options
You can control a variety of options for more advanced usage. Please look at the examples section for more information about how.
Here is a list of things that are available to control:
- use GPU or CPU.
- use color features.
- output the aligned 3D structures to SDF.
- direct input and output with RDkit molecules.
- pre-prepare Roshambo2 binary HDF5 files for fast runtime on large libraries.
- run an Roshambo2 server instance. 


## Color features
You can use color features for computing the score. These are pharmacophore features that are assigned to the molecules based on SMARTS matches. For each feature a dummy color atom is added to the 3D structure. These color features are then used to compute the score. To use the default colors you just need to pass the color=True flag to roshambo2 initialization and then use compute_method="combination" in the compute method:
```
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True)
scores = roshambo2_calculator.compute(optim_mode='combination') 

```
You can use custom color definitions by either:
- provide a custom RDKit compatible .fdef file.
- providing SMARTs for each feature.
- writing your own color feature generator (that can use 3rd party python code).
Please see the example: `example/custom_features.py`.

You can visualize the assigned color features to check they are what you expect. Please look at the example: `example/visualize_color.py`.



## Roshambo2 H5 files
The roshambo2 program can optionally take roshambo2 h5 files as input. These can be created using the preparation script. Reading data from these h5 files is much faster than reading from SDF files. The h5 file format contains groups of molecules. The group size is variable. Each group will be computed in one batch by the roshambo2 compute method. Within a molecule group the 3D coordinates of the molecules are stored as a h5 dataset that is a 3D array. The molecules are all padded to have the same number of atoms. The shape of the main coordinate array is [Number of configurations, number of atoms (padded), 4 ]. Where the 4th index is currently used to mark an atom as a padded atom. In future this could be used to store atom radii for the Gaussian overlap. Furthermore, there are dataset arrays for the molecule names, number of real atoms, number of color dummy atoms, and color type. Internally roshambo2 stores data in the same type of arrays.

If you have multiple roshambo2 h5 files you can combine them into one using the 'roshambo2.prepare.combine_datasets' function:
```python
from roshambo2.prepare import combine_datasets

datasets = ['processed_dataset_1.h5', 'processed_dataset_2.h5']

combine_datasets(datasets, 'combined.h5')
```
This will load the two datasets and combine them into one. It will also combine different molecule groups into one group. Note that you can set the approximate maximum size of the molecule groups. The default is 1 million.
```python
combine_datasets(datasets, 'combined.h5', max_mols_per_group=100000)
```
This can be useful if you have a GPU with small memory as currently the program will compute each molecule group as one batch on the GPU.

## Preparation script

The slowest part of the program is the assignment of color features using RDKit. The second slowest part is reading in the 3D SDF files. To ease the searching of large datasets we created a Roshambo2 data format (h5 file) and a script that will read in 3D SDF files, assign color features if requested, and create formatted Roshambo2 H5 files. The H5 files can be read in very quickly. The idea is you can prepare the dataset H5 files ahead of time and then run searches quickly each time you have a new query molecule.

```python
prepare_dataset_from_sdf.py --color dataset.sdf processed_dataset.h5
```

If the color flag is omitted then color features will not be assigned (much faster).


## Preparing Roshambo2 h5 files directly from RDKit molecules

It is possible to go directly from a list of RDKit molecules (with 3D conformations) and create an Roshambo2 h5 file.
```python
from roshambo2.prepare import prepare_from_rdkitmols

roshambo2_data = prepare_from_rdkitmols(dataset_mols, color=True)
    
  

roshambo2_data.save_to_h5('from_rdkit.h5')
```
For a full example look at "example/prepare_from_rdkit.py"


## Server mode
To setup server mode look at `examples/server_mode/`

You will need to create a config.yaml file by copying and modifying one of the examples.
e.g.:
```yaml
# config.yaml
dataset_files:
  - "../dataset.sdf"
hostname: "0.0.0.0"
port: 8087
api_name: "/roshambo2/search"
verbosity: 2

```

You can edit this script to load the datasets you want. You can also set the URL and port.
By default the IP is 0.0.0.0 which will run the server at `localhost`, and it will be accessible from other machines using the IP address of the host machine. (Providing the ports are open).

To run this example first make sure you are in the examples/server_mode folder (this is just to make sure the paths to the example dataset are correct)

Then run the roshambo2 server mode program.
```bash

roshambo2_server_app.py config.yaml

```
This python script should get installed to your path (the same as the prepare_from_sdf.py one). If for some reason is does not you can run run it
with python directly using its full path:
```bash
python <location of roshambo2 git repo>/roshambo2/scripts/roshambo2_server_app.py config.yaml
```

This will launch a Flask server that has loaded the roshambo2 dataset into memory. For a large dataset it will take some time to load it all in. It if it has loaded successfully the last output will look something like this:

```bash
...
INFO:roshambo2.roshambo2:Roshambo2 setup completed in 0.4868598530010786s
 * Serving Flask app 'server_app'
 * Debug mode: off
INFO:werkzeug:WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
 * Running on all addresses (0.0.0.0)
 * Running on http://127.0.0.1:8087
 * Running on http://192.168.1.69:8087
 ```

This means it has worked and it up and running.
This program will need to stay running.
You can leave it running in its own terminal, or run it in the background or similar.
For debugging it will be useful to see the output it prints to stdout and stderr.

Note that the warning message comes from the python Flask library. For internal network use you can ignore it. (If you wish to open the server to the web you will need to pay attention to it.)

Once this server program is running you can then send search queries to it. Please look at the example in `examples/server_mode/search_server.py`. Remember that the port and url used in that script will need to correspond to the port and url in the config.yaml.

It is recommended you first run the search_server.py script on the same machine that is running the server. Once you have checked that this is working correctly you can then open the required port to other machines in the network and run search queries from different machines.

 