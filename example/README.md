
# Examples

## 1. Basic run

look at `basic_run.py`. You can run it with:

```
python basic_run.py
```
the input is two sdf files: the query molecules and the dataset molecules

additionally to run with color features used in the optimizer look at:

`basic_run_color.py`


## 2. Preparing datasets before running and large scale runs

look at `large_scale_run.sh` first and then `large_run.py`
you can run with it with:
```
bash large_scale_run.sh
```
This shows how you prepare a dataset before.


## 3. Combine Roshambo2 h5 datasets

processed dataset files can be combined. Look at `combine_h5.py`

## 4. Server mode
Please look in the server_mode subfolder

```{include} ../example/server_mode/README.md
```


## 5. Modifying the optimization settings

`modify_optimizer_settings.py`


## 6. Integration with RDKit
`rdkit_integration.py` demonstrates how to use a list of RDKit molecules directly as the input.

`prepare_from_rdkit.py` demonstrates how to prepare an Roshambo2 H5 file from a list of RDKit molecules.

## 7. Customizing color/pharmacophore features
`custom_features.py`

## 8. visualizing the color features
`visualize_color.py`.
You can get the color atoms written as a pseudo-molecule that can be visualized with a standard molecule viewer, overlayed with the standard molecules:

![color_dummy_atoms](_img/color_dummy_atoms.png)
