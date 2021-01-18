# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hoenig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/usr/bin/env python3

from rdkit import Chem
import argparse


def convert(input, property, mol_id, output="result.smi"):
    sdf = Chem.SDMolSupplier(input)
    with open(output, 'w') as f:
        for mol in sdf:
            smi = Chem.MolToSmiles(mol)
            id = mol.GetProp(mol_id)
            prop = mol.GetProp(property)
            f.write(f"{smi}\t{id}\t{prop}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a SDF into a Smiles file.')
    parser.add_argument('input', metavar='INPUT', type=str, help='Input file in SD format.')
    parser.add_argument('property_col', metavar='PROPERTY_COLUMN', type=str, help='Column identifyer for molecule property of interest.')
    parser.add_argument('-i', '--id', metavar='ID_TAG', default='_Name', type=str, help="ID tag inside the input sdf file. If no tag is included, please use the default taking the molecules name of the header.")
    parser.add_argument('--output', metavar='output', type=str, help='Output file name')
    args = parser.parse_args()
    if args.output:
        convert(args.input, args.property_col, args.id, args.output)
    else:
        convert(args.input, args.property_col, args.id)
