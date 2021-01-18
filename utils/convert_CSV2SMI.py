# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/usr/bin/env python3

from os import path, system
import sys
import csv
import re
import argparse


# User: Specify CSV format below!
# key (string): property name, value (number): index of corresponding column inside the csv
# 'id' and 'smiles' are obligatory! Other property names can be user defined.
# If no id is available from CSV, specify it as 'nan' to get an id from a counter.
format = {'id': 1, 'mol_weight': 5, 'alogp': 6, 'smiles': 0, 'pic50': 4}


# Example formats:
# ChEMBL CSV: {'id': 0, 'mol_weight': 3, 'alogp': 5, 'smiles': 7, 'ic50': 10}
# MoleculeNet data sets:
# BACE CSV: {'id': 1, 'mol_weight': 5, 'alogp': 6, 'smiles': 0, 'pic50': 4}
# ESOL CSV: {'id': "nan", 'mol_weight': 3, 'esol': 1, 'smiles': 9, 'sol': 8, 'surface': 7}
# FreeSolv: {'id': "nan", 'smiles': 1, 'expt': 2, 'calc': 3}
# Lipophilicity: {'id': 0, 'smiles': 2, 'logD': 1}

def convert_csv2smiles(file_path, delim, out_name="result.smi"):
    """
    Converts a CSV file in above specified format into a SMILES file fitting the format of mmpdb.
    Saves all important information in an additional dictionary named molecules.
    Keys of molecules are the IDs.
    @param file_path: The CSV input file path.
    @param delim: The CSV delimiter.
    @param out_name: The name of the output smiles file to be written.
    @return: The path to the output smiles file and the molecules array.
    """
    # Molecule is identified by its smiles
    molecules = {}
    try:
        with open(file_path, newline='') as csvfile:
            mol_reader = csv.reader(csvfile, delimiter=delim, quotechar='"')
            next(mol_reader)  # Remove row containing identifiers
            for i, row in enumerate(mol_reader):
                # Specify columns of id, smiles and other properties you like to analyze inside the CSV by adjusting the properties dict.
                # Adjust numbers inside row[*] to fit to columns inside your CSV and name the properties of interest as below.
                # id and smiles are obligatory! If no id is available, replace row[*] by i for an id counter.
                properties = {}
                for key in format:
                    if format[key] == 'nan':
                        properties[key] = f"{i}"
                    else:
                        properties[key] = row[format[key]]
                molecules[properties['id']] = properties
        output = open(out_name, 'w')
        for id in molecules:
            output.write(f"{molecules[id]['smiles']}\t{id}\t{molecules[id]['pic50']} \n")
        output.close()
        print(f"Your file was parsed into  {out_name} (Smiles format)")
        return out_name, molecules
    except FileNotFoundError:
        print(f"Your file was not found. No such file or directory: {file_path}")
        print("Please check your input and try again.")
        sys.exit(1)
    except:
        print("Your file has the wrong format.")
        print("Please adapt the properties dict to fit to your CSV layout.")
        sys.exit(1)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a SDF into a Smiles file.')
    parser.add_argument('input', metavar='INPUT', type=str, help='Input file in SD format.')
    parser.add_argument('property', metavar='PROPERTY', type=str, help='Property of interest (e.g. pic50).')
    parser.add_argument('--delimiter', metavar='DELIMITER', default=",", type=str, help='CSV delimiter. Default: ",".')
    parser.add_argument('-o', '--output', metavar='OUTPUT', type=str, help='Output file name')
    args = parser.parse_args()
    if args.output:
        convert_csv2smiles(args.input, args.delimiter, args.output)
    else:
        convert_csv2smiles(args.input, args.delimiter)