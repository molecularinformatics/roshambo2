
from roshambo2 import Roshambo2

# we can customize the features by using the PharmacophoreGenerator


from roshambo2.pharmacophore import RDKitPharmacophoreGenerator, CustomSMARTSPharmacophoreGenerator


## Example using default RDKit features

# 1.  create the list of features, this example each features uses the RDKit definition
features = {'Donor': 'rdkit',
        'Acceptor': 'rdkit',
        'PosIonizable': 'rdkit',
        'NegIonizable': 'rdkit',
        'Aromatic': 'rdkit',
        'Hydrophobe': 'rdkit',
    }  # these are also the default features

# the available default features are those found in https://github.com/rdkit/rdkit/blob/master/Data/BaseFeatures.fdef
# you can provide your own fdef file as `fdefname='your_file.fdef'` when you create the PharmacophoreGenerator class below. (Not yet tested)


# 2. create the interaction matrix
# the interaction map has an entry for each feature-feature interaction you want. 
# The 3rd and 4th entries are the gaussian width and height respectively
# If an interaction is not in the map then it will be 0,0 i.e. no interaction.
interactions = [
            ('Donor', 'Donor', 1.0, 1.0),
            ('Acceptor', 'Acceptor', 1.0, 1.0),
            ('PosIonizable', 'PosIonizable', 1.0, 1.0),
            ('NegIonizable', 'NegIonizable', 1.0, 1.0),
            ('Aromatic', 'Aromatic', 1.0, 1.0),
            ('Hydrophobe', 'Hydrophobe', 1.0, 1.0)
        ] # there are also the default interactions

# 3. create the color generator
color_generator =  RDKitPharmacophoreGenerator(features=features, interactions=interactions)
# we can check the indexes that correspond to each feature:
print(color_generator.get_feature_indexes())

roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True, verbosity=1, color_generator = color_generator)

scores = roshambo2_calculator.compute(backend='cpp', optim_mode='combination', reduce_over_conformers=False, start_mode=1) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)



## Example using custom RDKit fdef file
# we can use reuse the features and interaction map from before, but this time we use a specific RDKit fdef file
color_generator =  RDKitPharmacophoreGenerator(features=features, interactions=interactions, fdefName='./minimal.fdef')

roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True, verbosity=1, color_generator = color_generator)

scores = roshambo2_calculator.compute(backend='cpp', optim_mode='combination', reduce_over_conformers=False, start_mode=1) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)


## Example using custom SMARTS (via RDKit)

# To use custom smarts to be matched with RDkit the features dict needs to have the values as a list of the smarts patterns, e.g.:
features= {
  "Donor": ["[#16!H0]"],
  "Acceptor": ["[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"],
  "PosIonizable": ["[+,+2,+3,+4]"],
  "NegIonizable": ["[-,-2,-3,-4]"],
  "Aromatic": ["a1aaaaa1", "a1aaaa1"],
  "Hydrophobe": ["[C&r3]1~[C&r3]~[C&r3]1"],   
}

# you will need to supply the interaction matrix also
interactions = [
            ('Donor', 'Donor', 1.0, 1.0),
            ('Acceptor', 'Acceptor', 1.0, 1.0),
            ('PosIonizable', 'PosIonizable', 1.0, 1.0),
            ('NegIonizable', 'NegIonizable', 1.0, 1.0),
            ('Aromatic', 'Aromatic', 1.0, 1.0),
            ('Hydrophobe', 'Hydrophobe', 1.0, 1.0)
        ] 

color_generator = CustomSMARTSPharmacophoreGenerator(features=features, interactions=interactions)

roshambo2_calculator = Roshambo2("query.sdf", "dataset.sdf", color=True, verbosity=1, color_generator = color_generator)

scores = roshambo2_calculator.compute(backend='cpp', optim_mode='combination', reduce_over_conformers=False, start_mode=1) 

for key in scores:
    print(key)
    df = scores[key]
    print(df)



## Full customization
# To use a different program to match the smarts and generate the coordinate of the dummy atoms you can create your own PharmacophoreGenerator class
# Use the roshambo2/pharmacophore.py CustomSMARTSPharmacophoreGenerator class as your example.