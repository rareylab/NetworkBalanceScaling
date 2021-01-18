# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/usr/bin/env python3

from rdkit import Chem
import csv
import numpy as np

format = {'id': 1, 'mol_weight': 5, 'alogp': 6, 'smiles': 0, 'pic50': 4}
delimiter = ","
input = "bace.csv"
train = "bace_train.sdf"
test = "bace_test_with_properties.sdf"
prop_of_interest = "pic50"
# Example formats:
# ChEMBL CSV: {'id': 0, 'mol_weight': 3, 'alogp': 5, 'smiles': 7, 'ic50': 10}
# MoleculeNet data sets:
# BACE CSV: {'id': 1, 'mol_weight': 5, 'alogp': 6, 'smiles': 0, 'pic50': 4}
# ESOL CSV: {'id': "nan", 'mol_weight': 3, 'esol': 1, 'smiles': 9, 'sol': 8, 'surface': 7}
# FreeSolv: {'id': "nan", 'name': 0,'smiles': 1, 'expt': 2, 'calc': 3}
# Lipophilicity: {'id': 0, 'smiles': 2, 'logD': 1}
mols = []

with open(input, newline='') as csvfile:
    mol_reader = csv.reader(csvfile, delimiter=delimiter, quotechar='"')
    next(mol_reader)  # Remove row containing identifiers
    for info in mol_reader:
        m = Chem.MolFromSmiles(info[format["smiles"]])
        for num, key in enumerate(format.keys()):
            if format[key] == "nan":
                m.SetProp(key, str(num))
            else:
                m.SetProp(key, info[format[key]])
        mols.append(m)
train_out = Chem.SDWriter(train)
test_out = Chem.SDWriter(test)

# Prepare data split
folds = 10
np.random.seed(1234)
np.random.shuffle(mols)
subsets = np.array_split(mols, folds)

# For a k-fold validation, create a loop to consider all subsets
# Here, only the first subset is split as a test set
for mol in mols:
    if mol in subsets[0]:
        # Reactivate line below to delete properties of test set molecules before writing them
        #mol.SetProp(prop_of_interest, "nan")
        test_out.write(mol)
    else:
        train_out.write(mol)
