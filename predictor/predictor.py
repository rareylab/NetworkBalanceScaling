# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey
 
# Parts of the DeepChem library (https://github.com/deepchem/) are included.
# Thus, this script is additionally licensed unter MIT License.

#!/usr/bin/env python3

from rdkit import Chem
import deepchem
import argparse
from deepchem.utils.save import load_dataset_from_disk
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="predictor.py QUERY [-d MODEL_DIR, -o RESULT_NAME]")
    parser.add_argument("query", metavar="QUERY", type=str, help="Query molecules in SD format.")
    parser.add_argument("task", metavar="TASK", type=str,
                        help="Task for the model to perform. Depends on the input data set.",
                        choices=["pic50", "logD", "expt", "sol"])
    parser.add_argument("-d", "--modeldir", metavar="MODEL_DIR", default="model_dir", type=str,
                        help="Location, where the predictor was stored.")
    parser.add_argument("-o", "--outfile", metavar="RESULT_NAME", default="prediction.smi", type=str,
                        help="Name of the result file in smi format.")
    parser.add_argument("-m", "--modeltype", metavar="MODEL_TYPE", default="graphconvreg",
                        choices=["graphconvreg", "mpnn"])
    parser.add_argument("-i", "--id", metavar="ID_TAG", default="_Name", help="ID tag inside the input sdf file. If no tag is included, please use the default taking the molecules name of the header.")
    args = parser.parse_args()

    tasks = [args.task]
    query_molecules = args.query
    model_type = args.modeltype

    # Restore model and featurizers
    if model_type == "graphconvreg":
        model = deepchem.models.GraphConvModel(
            len(tasks),
            model_dir=args.modeldir,
            mode='regression')
        featurizer_func = deepchem.feat.ConvMolFeaturizer()

    elif model_type == "mpnn":
        model = deepchem.models.MPNNModel(
            len(tasks),
            n_atom_feat=75,
            n_pair_feat=14,
            n_hidden=75,
            model_dir=args.modeldir,
            mode='regression')
        featurizer_func = deepchem.feat.WeaveFeaturizer()

    else:
        print('Unknown model type. Please select from ["graphconvreg", "mpnn"]')
        sys.exit(1)
    model.restore()

    # Load data set and transformers
    loader = deepchem.data.SDFLoader(tasks=tasks, clean_mols=True, featurizer=featurizer_func)
    query_dataset = loader.featurize(query_molecules)
    loaded, model_datasets, transformers = load_dataset_from_disk(args.modeldir)
    if not loaded:
        print("Error loading transformers. Please try again.")
        sys.exit(1)
    for transformer in transformers:
        transformer.transform(query_dataset)

    # Predict
    predictions = [x[0] for x in model.predict(query_dataset, transformers)]

    # Save predictions
    ids = query_dataset.ids
    zipped = zip(ids, predictions)
    dictionary = dict(zipped)
    keys = dictionary.keys()
    result = open(args.outfile, "w")
    suppl = Chem.SDMolSupplier(query_molecules)
    for mol in suppl:
        smiles = Chem.MolToSmiles(mol)
        id = mol.GetProp(args.id)
        if smiles in keys:
            result.write(f"{smiles}\t{id}\t{dictionary[smiles]}\n")
        else:
            print(f"WARNING: Molecule {id} was not predicted!\n")
    result.close()
    print(f"Predictions were written to {args.outfile} in format SMILES\\tID\\tPROPERTY\\n.")
