
from roshambo2 import Roshambo2
 
roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True, verbosity=2, remove_Hs_before_color_assignment=False)

scores = roshambo2_calculator.compute(backend='cpp', optim_mode='combination', combination_param=0.5, reduce_over_conformers=False, start_mode=1, keep_order=False) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)



# setting write_color_pseudomols=True means a second SDF file will be written. 
# this contains pseudmolecules corresponding to the molecules in the main SDF file.
# The atoms are the feature color dummy atoms. Their symbols use this default mapping:
#{'feature1':'X1', 'feature2':'X2', ... ,'featureN': 'XN'}

#roshambo2_calculator.write_best_fit_structures(write_color_pseudomols=True, append_query=True)


# You can assign what symbols you want (e.g. real elements) to color features so we can easily see them in a standard molecular viewer
# you provide a dictionary that maps from feature to symbol, e.g.
feature_to_symbol = {'Donor':'H', 'Acceptor':'He', 'PosIonizable':'Li', 'NegIonizable':'Be', 'Aromatic':'B', 'Hydrophobe':'C'}
roshambo2_calculator.write_best_fit_structures(write_color_pseudomols=True, append_query=True, feature_to_symbol_map=feature_to_symbol)