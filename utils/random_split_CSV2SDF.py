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
train = "bace_train_tvt.sdf"
test = "bace_test_tvt.sdf"
validation = "bace_validation_tvt.sdf"
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
validation_out = Chem.SDWriter(validation)

# Prepare data split
folds = 10
np.random.seed(1234)
np.random.shuffle(mols)
train_set, validation_set, test_set = np.split(mols, [int(0.8 * len(mols)), int(0.9 * len(mols))]) # Generates an 80 - 10 - 10 split
for mol in mols:
    if mol in train_set:
        # Deactivate line below to not delete properties of test set molecules before writing them
        mol.SetProp(prop_of_interest, "nan")
        train_out.write(mol)
    elif mol in validation_set:
        validation_out.write(mol)
    else:
        test_out.write(mol)
print(len(train_set), len(validation_set), len(test_set))
# BACE: 1210 - 151 - 152
