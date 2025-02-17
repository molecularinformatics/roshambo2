# MIT License
# 
# Copyright (c) 2025 molecularinformatics  
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import os

import numpy as np

from rdkit.Chem import AllChem, MolFromSmarts
from rdkit.Chem.rdchem import Mol as rdkitMol
from rdkit import RDConfig

class BasePharmacophoreGenerator():
    """Base Class for Roshambo2 PharmacophoreGenerators

       Functional PharmacophoreGenerators classes are subclass of this.

       If you want to make your own please see the user/developer guide.

    """
    def __init__(self):
        self.FEATURES = {}

        # enumerate for better compatibility with backend codes
        # Note count from 1, type 0 are normal atoms
        self.FEATURES_ENUM = { key: i+1 for i,key in enumerate(self.FEATURES.keys())}



    def set_interactions(self, FEATURES_INTERACTIONs):
        """set the interaction matrix.

            Args: 
                FEATURES_INTERACTIONs (List[Tuple]): map of the interactions between pairs of features, 
                    e.g. [ ('feature1', 'feature2', 1.0, 1.0), ... ]. Where the last two values are the gaussian width and height. 
                    Any pair interactions not included are set to zero. 
        
        """

        # TODO: checks
        self.FEATURES_INTERACTIONs = FEATURES_INTERACTIONs
        
        # turn it into a matrix, note zeros for non interactions
        self.interaction_matrix_r = np.ones((len(self.FEATURES)+1, len(self.FEATURES)+1))
        self.interaction_matrix_p = np.zeros((len(self.FEATURES)+1, len(self.FEATURES)+1))

        for interaction in self.FEATURES_INTERACTIONs:
            key1 = interaction[0]
            key2 = interaction[1]
            r = interaction[2]
            p = interaction[3]

            i = self.FEATURES_ENUM[key1]
            j = self.FEATURES_ENUM[key2]

            self.interaction_matrix_r[i,j] = r
            self.interaction_matrix_p[i,j] = p
        
        # print(self.interaction_matrix_r)
        # print(self.interaction_matrix_p)

    def get_feature_indexes(self):
        """Get the mapping from feature name to index.
            
            Returns:
                Dict: Mapping from feature name (str) to index (int) counting starts at 1.
        """

        return self.FEATURES_ENUM
    
    def get_index_to_feature(self):

        return {i:F for F,i in self.FEATURES_ENUM.items()}

    def generate_color_atoms(self, rdkitMol):
        """ takes an RDKit molecule and returns the color atoms

            you will need to have created an assign_features method for you subclass.
        
            return type should be a list of coordinates and list of the color type
            the color type is an integer corresponding to the self.FEATURES_ENUM of the class.

            Args:
                rdkitMol (rdkitMol): An RDKit molecule.

            Returns:
                np.array, dtype=float32: Numpy array of shape [l,3] containing the coordinates of l color dummy atoms.
                np.array, dtype=int: Numpy array of shape [l] containing the color type of each color dummy atom. The indexes must correspond to self.FEATURES_ENUM.
        """
    
        features = self.assign_features(rdkitMol)
        l=len(features)
        color_coords = np.zeros((l,3),dtype=np.float32)
        color_type = np.zeros(l, dtype=int)

        for i in range(l):
            feature = features[i]
            color_coords[i,:] = np.array(feature[2], dtype=np.float32)
            color_type[i] = self.FEATURES_ENUM[feature[0]]


        return color_coords, color_type





class RDKitPharmacophoreGenerator(BasePharmacophoreGenerator):
    """Assigns pharmacophore features to an RDKit molecule using RDKit default fdef files.

    This class is modified from Roshambo provided under MIT license:
    https://github.com/molecularinformatics/roshambo_biogen/tree/Paperless

    """
    def __init__(self, features = None, interactions = None, fdefName = None):
        """
            Args:
                features (Dict, Optional): Feature dictionary, if None the default is used.
                interactions (List, Optional): Interaction map, if None the default is used.
                fdefName (str, Optional): RDKit fdef file to use, if None the default RDKit BaseFeatures.fdef is used.

        """


        if fdefName is None:
            self.fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        else:
            self.fdefName = fdefName

        # self.factory = AllChem.BuildFeatureFactory(self.fdefName)
        # wait untill we need it to build it because it cannot be pickled by python
        self.factory = None

        if features is None:
            self.FEATURES = {
                'Donor': 'rdkit',
                'Acceptor': 'rdkit',
                'PosIonizable': 'rdkit',
                'NegIonizable': 'rdkit',
                'Aromatic': 'rdkit',
                'Hydrophobe': 'rdkit',
            }
        else:
            self.FEATURES = features


        # enumerate for better compatability with backend codes
        # Note count from 1, type 0 are normal atoms
        self.FEATURES_ENUM = { key: i+1 for i,key in enumerate(self.FEATURES.keys())}

        # define the interactions
        # format (type, type, radius, height)
        if interactions == None:
            FEATURES_INTERACTIONs = [
                ('Donor', 'Donor', 1.0, 1.0),
                ('Acceptor', 'Acceptor', 1.0, 1.0),
                ('PosIonizable', 'PosIonizable', 1.0, 1.0),
                ('NegIonizable', 'NegIonizable', 1.0, 1.0),
                ('Aromatic', 'Aromatic', 1.0, 1.0),
                ('Hydrophobe', 'Hydrophobe', 1.0, 1.0)
            ]
        else:
            FEATURES_INTERACTIONs = interactions

        self.set_interactions(FEATURES_INTERACTIONs)


    def assign_features(self, mol: rdkitMol):
        """Assigns features to an RDKit molecule.
        
            Args:
                mol (rdkitMol)

            Returns:
                List: A list of the features. Will be used by the generate_color_atoms method.
        
        """
        # first time:
        if self.factory is None:
            # print("creating RDKit FeatureFactory")
            self.factory = AllChem.BuildFeatureFactory(self.fdefName)
        
        features = self.factory.GetFeaturesForMol(mol)

        pharmacophores = []
        for feature in features:
            fam = feature.GetFamily()
            if fam in self.FEATURES.keys():
                pos = feature.GetPos()
                atom_indices = feature.GetAtomIds()
                p = [
                    fam,
                    atom_indices,
                    [pos[0], pos[1], pos[2]],
                ]
                pharmacophores.append(p)
        return pharmacophores




class CustomSMARTSPharmacophoreGenerator(BasePharmacophoreGenerator):
    """Assigns pharmacophore features to an RDKit molecule using custom SMARTS via RDKit.

    This class is modified from Roshambo provided under MIT license:
    https://github.com/molecularinformatics/roshambo_biogen/tree/Paperless

    """
    def __init__(self, features = None, interactions = None):
        """
            Example:
                see /example/customs_features.py
                
            Args:
                features (Dict): Dictionary where the key is the feature name and the value is a SMARTS pattern that defines the feature.
                interactions (List): An interaction map for feature pair interactions. Any pair not specified will have their interactions set to zero.
        """

        # no defaults here, user must supply both features and interactions
        assert(features is not None and interactions is not None)
        self.FEATURES = features
        self.FEATURES_ENUM = { key: i+1 for i,key in enumerate(self.FEATURES.keys())}

        self.compiled_smarts = {k: list(map(MolFromSmarts, v)) for k, v in self.FEATURES.items()}

        self.set_interactions(interactions)



    def compute_match_centroid(self, mol, matched_pattern):
        conf = mol.GetConformer()
        positions = [conf.GetAtomPosition(i) for i in matched_pattern]
        center = np.mean(positions, axis=0)
        return tuple(center)


    def find_matches(self, mol, patterns):
        matches = []
        for pattern in patterns:
            # Get all matches for that pattern
            matched = mol.GetSubstructMatches(pattern)
            for m in matched:
                # Get the centroid of each matched group
                # centroid = average_match(mol, m)
                centroid = self.compute_match_centroid(mol, m)
                # Add the atom indices and (x, y, z) coordinates to the list of matches
                matches.append([m, centroid])
        return matches


    def assign_features(self, rdkit_mol):
        matches = {}
        for key, value in self.compiled_smarts.items():
            matches[key] = self.find_matches(rdkit_mol, value)

        # Sometimes, a site can match multiple SMARTS representing the same pharmacophore,
        # so we need to keep it only once
        cleaned_matches = {}
        for key, value in matches.items():
            unique_lists = []
            for lst in value:
                if lst not in unique_lists:
                    unique_lists.append(lst)
            cleaned_matches[key] = unique_lists

        pharmacophore = []
        for key, value in cleaned_matches.items():
            for match in value:
                p = [key, match[0], match[1]]
                pharmacophore.append(p)

        return pharmacophore
        


# Default
PharmacophoreGenerator=RDKitPharmacophoreGenerator
