#!/usr/bin/env python

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

import argparse
import logging
from roshambo2.prepare import prepare_from_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create Roshambo2 h5 dataset")
    parser.add_argument("input_file", help="input sdf file containing molecules with 3D coordinates")
    parser.add_argument("output_file", default="roshambo2_dataset.h5", help="Output file name. Default is 'roshambo2_dataset.h5.h5'")
    parser.add_argument("--conformers_have_unique_names", default=False, type=bool, help="If true it assumes you have named the conformers of MOLA as MOLA_0, MOLA_1 etc. If false it assumes conformers of MOLA are all called MOLA and the program will append the _N to distinguish them.")
    parser.add_argument("-c", "--color", action="store_true", help="Setup color features")
    parser.add_argument("-m", "--max_mols_per_group", default=1000000, help="max number of molecules per group in Roshambo2 H5 dataset." )
    parser.add_argument("-v", "--verbosity", default=1, type=int, help="verbosity level options are [0,1,2].")
    parser.add_argument("-n", "--n_cpus", default=None, type=int, help="number of cpus to use (default is all)")
    args = parser.parse_args()
    print(args)
    # TODO enable custom color features here

    if args.verbosity==0:
        quiet=True # dont show the progress bars
    else:
        quiet=False


    assert(args.verbosity in [0,1,2])

    logging.basicConfig(
        level=[logging.WARNING, logging.INFO, logging.DEBUG][args.verbosity],
        # format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        # datefmt='%Y-%m-%d %H:%M:%S'
    )

    
    prepare_from_file(args.input_file, args.output_file, color=args.color, conformers_have_unique_names=args.conformers_have_unique_names, max_mols_per_group=int(args.max_mols_per_group), quiet=quiet, n_cpus=args.n_cpus)