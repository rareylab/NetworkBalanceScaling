#!/usr/bin/env python3

# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey



from rdkit import Chem
import argparse
import math
import os
from sklearn.metrics import mean_squared_error, r2_score
from statistics import stdev

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the RMSD value for a molecular property.')
    parser.add_argument('--original', '-o', metavar='ORIGINAL', required=True, type=str, help='File with original property values. SDF or SMILES format.')
    parser.add_argument('--prediction', '-p', metavar='PREDICTION', required=True, type=str, help='File with predicted property values. SDF or SMILES format.')
    parser.add_argument('--optimized', '-q', metavar="OPTIMIZED", type=str, help="File with optimized property predictions. SDF or SMILES format.")
    parser.add_argument('-i', '--id', type=str, default="_Name", help="ID tag for files in SD format.")
    parser.add_argument('-l', '--logfile', type=str, default="evaluation_result.log")
    parser.add_argument('--label_predicted', type=str, required=True, help="Label of the prediction.")
    parser.add_argument('--label_optimized', type=str, required=True, help="Label of the optimized prediction.")
    
    parser.add_argument('property', metavar='Property', type=str, help='Property of interest (should be similar in both files.')
    args = parser.parse_args()

    original_dict = {}
    pred_dict = {}
    opti_dict = {}
    orig_ext = os.path.splitext(args.original)[-1].lower()
    pred_ext = os.path.splitext(args.prediction)[-1].lower()

    id_label = args.id

    # collect original values
    if orig_ext == ".sdf":
        original_mols = Chem.SDMolSupplier(args.original)
        for mol in original_mols:
            original_dict[mol.GetProp(id_label)] = float(mol.GetProp(args.property))
    else:
        for line in open(args.original, "r"):
            line_arr = line.split("\t")
            original_dict[line_arr[1]] = float(line_arr[2])

    # collect values predicted values
    no_deviation = False
    if pred_ext == ".sdf":
        pred_mols = Chem.SDMolSupplier(args.prediction)
        for mol in pred_mols:
            pred_dict[mol.GetProp(id_label)] = float(mol.GetProp(args.property))
    else:
        for line in open(args.prediction, "r"):
            line_arr = line.split("\t")
            pred_dict[line_arr[1]] = float(line_arr[2])

    if args.optimized:
        if os.path.splitext(args.optimized)[-1].lower() == ".sdf":
            pred_mols = Chem.SDMolSupplier(args.optimized)
            for mol in pred_mols:
                opti_dict[mol.GetProp(id_label)] = float(mol.GetProp(args.label_optimized))
        else:
            for line in open(args.prediction, "r"):
                line_arr = line.split("\t")
                pred_dict[line_arr[1]] = float(line_arr[2]) + float(line_arr[3])

    sum_sq = 0

    preds = []
    orgs = []
    optis_all = [] # including unoptimizable prediction values
    optis_only = []
    orgs_only = []
    unoptimizable_ids = []
    
    unopt_values = []
    unopt_values_for_ori = []

    pred_in_net = []
    pred_not_in_net = []

    for id in pred_dict.keys():
        preds.append(pred_dict[id])
        orgs.append(original_dict[id])
        diff = (pred_dict[id] - original_dict[id])
        if args.optimized:
            if id in opti_dict.keys():
                optis_only.append(opti_dict[id])
                orgs_only.append(original_dict[id])
                optis_all.append(opti_dict[id])
                pred_in_net.append(pred_dict[id])
            else:
                unoptimizable_ids.append(id)
                optis_all.append(pred_dict[id])

                unopt_values.append(pred_dict[id])
                unopt_values_for_ori.append(original_dict[id])

                pred_not_in_net.append(pred_dict[id])

    with open(args.logfile, "w") as f:
        stdDev_orig_all = stdev(orgs)
        stdDev_orig_only_in_net = stdev(orgs_only)
        stdDev_orig_only_not_in_net = stdev(unopt_values_for_ori)

        stdDev_opt_all = stdev(optis_all)
        stdDev_opt_only_in_net = stdev(optis_only)

        f.write("StdDevs:\n")
        f.write(f"original all values:              {stdDev_orig_all}\n")
        f.write(f"original only in net values:      {stdDev_orig_only_in_net}\n")
        f.write(f"original only NOT in net values:  {stdDev_orig_only_not_in_net}\n")

        f.write(f"optimized all values:             {stdDev_opt_all}\n")
        f.write(f"optimized only in net values:     {stdDev_opt_only_in_net}\n\n")


        f.write(f"                              Root Mean Square Deviation    R²        Num_of_Samples\n")
        f.write(f"Predictions(all):             {mean_squared_error(orgs, preds, squared=False):.6f}                      {r2_score(orgs, preds):.6f}  {len(orgs)}\n")
        f.write(f"Predictions(only in net):     {mean_squared_error(orgs_only, pred_in_net, squared=False):.6f}                      {r2_score(orgs_only, pred_in_net):.6f}  {len(orgs_only)}\n")
        f.write(f"Predictions(only NOT in net): {mean_squared_error(unopt_values_for_ori, pred_not_in_net, squared=False):.6f}                      {r2_score(unopt_values_for_ori, pred_not_in_net):.6f}  {len(unopt_values_for_ori)}\n\n")

        if len(optis_all) > 0:
            f.write(f"Optimized (all):              {mean_squared_error(orgs, optis_all, squared=False):.6f}                      {r2_score(orgs, optis_all):.6f}  {len(orgs)}\n")
            f.write(
                f"Optimized (only):             {mean_squared_error(orgs_only, optis_only, squared=False):.6f}                      {r2_score(orgs_only, optis_only):.6f}  {len(orgs_only)}\n\n")
        if len(unopt_values) > 0:
            f.write(f"Scores (un_opt):              {mean_squared_error(unopt_values_for_ori, unopt_values, squared=False):.6f}                      {r2_score(unopt_values_for_ori, unopt_values):.6f}  {len(unopt_values)}\n")
        f.write(f"\nUnoptimizable molecules (IDs) ({len(unoptimizable_ids)} mols):\n")
        for id in unoptimizable_ids:
            f.write(f"{id}\n")
