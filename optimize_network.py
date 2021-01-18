# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/usr/bin/env python3

import argparse
from os import path, system
import sys
import igraph as ig
import cvxpy as cp  # For convex optimization
import math
from rdkit import Chem


def optimize_network(g, num_trafo, std_dev, input_name, outdir, loga, prediction=False):
    """
    Calculate the deviation on the edges and vertices with unknown or already predicted property.
    This uses cvxpy as a convex solver to optimize the overall (sum of squares) deviation of the network.
    @param g: The network (igraph) to be optimized
    @param num_trafo: Indicates how strongly the number of occurrences of a transformation is included into the
                      corresponding edge's weight in the objective function of the network optimization
    @param std_dev: Indicates how strongly the standard deviation of the transformation at an edge is included into the
                    edge's weight in the objective function of the network optimization
    @param input_name: Used for internal file naming
    @param outdir: Directory name for resulting graph file and log
    @param loga: Decide whether to use logarithmic scaling for the number of transformation occurrences in the weighting
    @param prediction: Decide whether to optimize existing predictions or predict unknown properties
    """
    sum_sq = 0
    num_edges = len(g.es)
    num_red_nodes = len(g.vs.select(color='red'))
    total_variables = num_red_nodes + num_edges
    # Initialize space for variables of nodes without defined property
    g.vs.select(color="red")["variable"] = [-1] * num_red_nodes
    x = cp.Variable(total_variables)  # One variable for each edge and each unlabeled node
    y = cp.Variable(total_variables)  # This resembles A@x in CVXPY's syntax, i.e. The matrix to transform the equations
    constraints = []  # This will hold all constructed constraint equations of the system
    var_num = 0  # Provides the next free solver variable
    print('Constructing Constraints for simplex solver. This may take a few minutes...')
    blacklist = []
    for edge in g.es:
        # early exit, if edge's antiparallel sibling is already inside the network and blacklisting
        if edge in blacklist:
            continue
        source, target = edge.vertex_tuple
        blacklist.append(g.es.select(_source=target.index, _target=source.index))

        edge["variable"] = var_num
        if std_dev > 0 or num_trafo > 0:
            if loga:
                factor = std_dev * (1 / (1 + edge['std_dev'])) + num_trafo * math.log(1 + edge['transform_num'])
            else:
                factor = std_dev * (1 / (1 + edge['std_dev'])) + num_trafo * edge['transform_num']
        else:
            factor = 1
        constraints.append(y[var_num] == factor * x[var_num])
        var_num += 1
        # Assign variables to the vertices, if they do not already have them or do not have a property
        if source["color"] == 'red' and source["variable"] < 0:
            source["variable"] = var_num
            constraints.append(y[var_num] == 0 * x[var_num])
            var_num += 1
        if prediction:
            source_var = (x[source["variable"]] + source['property']) if source["color"] == 'red' else source["property"]
        elif source["color"] == 'red' and source["property"] != 0:
            print(
                "WARNING: Identified red vertex with property != 0.\nPlease choose option [-p] for networks containing existing predictions at red vertices.")
        else:
            source_var = x[source["variable"]] if source["color"] == 'red' else source["property"]
        if target["color"] == 'red' and target["variable"] < 0:
            target["variable"] = var_num
            constraints.append(y[var_num] == 0 * x[var_num])
            var_num += 1
        if prediction:
            target_var = (x[target["variable"]] + target['property']) if target["color"] == 'red' else target["property"]
        elif target["color"] == 'red' and target["property"] != 0:
            print(
                "WARNING: Identified red vertex with property != 0.\nPlease choose option [-p] for networks containing existing predictions at red vertices.")
        else:
            target_var = x[target["variable"]] if target["color"] == 'red' else target["property"]
        constraints.append(target_var - source_var == edge["transform_mean"] + x[edge["variable"]])
    objective = cp.Minimize(cp.sum_squares(y))

    prob = cp.Problem(objective, constraints)
    prob.solve()

    # Adjust edge labels
    for edge in g.es:
        # Overwrite the LP variable id with its value inside the igraph
        edge["variable"] = x.value[edge["variable"]]
        sum_sq += edge['variable'] * edge['variable']
        edge["label"] = f"{edge['transform_mean']:.3f}\ndeviation: {edge['variable']:.3f}"

    # Adjust vertex labels
    for vertex in g.vs:
        if vertex["color"] == 'red':
            if prediction:
                vertex["variable"] = x.value[vertex["variable"]]
                vertex["label"] = f"{vertex['property']},\ndeviation: {vertex['variable']}"
            else:
                vertex["property"] = x.value[vertex["variable"]]
                vertex["label"] = vertex["label"].replace("?", f"{vertex['property']}")

    # Save results inside another GraphML file
    if not path.isdir(f'{outdir}'):
        system(f'mkdir {outdir}')

    graph_name = path.join(outdir, f"{input_name}_optimized.GraphML")
    g.write_graphml(graph_name)
    log = open(f"{path.join(outdir, input_name)}_optimization_log.txt", 'w')
    log.write(f"File: {graph_name}, parameters: t = {num_trafo}, s = {std_dev}, logarithmic = {loga}\n")
    log.write(f"Total sum of squares deviation: {prob.value:.5f}\n")
    log.write(f"Total sum of squares without factors: {sum_sq:.5f}\n")
    log.write(f"Sum of Squares deviation per edge: {prob.value / len(g.es):.5f}\n")
    log.write(f"Sum of Squares deviation per edge without factors: {(sum_sq / len(g.es)):.5f}\n")
    log.write(f"Number of edges: {len(g.es)}\n")
    log.close()
    print(f"Optimization details can be found in {path.join(outdir, input_name)}_optimization_log.txt.")
    print(f"Your optimized network was saved in GraphML format in {graph_name}.")


def sdf_to_dict(sdf, id_label, prop_label):
    """
    Reads molecules and returns a dict[label] -> property-value.
    """
    result_dict = dict()
    for mol in Chem.SDMolSupplier(sdf):
        result_dict[mol.GetProp(id_label)] = mol.GetProp(prop_label)
    return result_dict
    

def write_results_smi(g, input_name, outdir, prediction=False):
    """
    Writes the predictions made during optimization into a smiles file.
    @param g: igraph to be evaluated.
    @param input_name: Name of the input file to adapt the name of the written output.
    @param outdir: Name of the output directory.
    @param prediction: Decide whether there are vertex deviations annotated or not.
    """
    results = open(path.join(outdir, f"{input_name}.smi"), "w")
    for vertex in g.vs.select(color='red'):
        try:
            if prediction:
                results.write(f"{vertex['name']}\t{vertex['chem_id']}\t{vertex['property']}\t{vertex['variable']}\n")
            else:
                results.write(f"{vertex['name']}\t{vertex['chem_id']}\t{vertex['property']}\n")
        except:
            sys.exit("Error writing results. Please check your input and try again.")
    results.close()
    if prediction:
        print(
            f"Gathered predictions were saved in {outdir}/{input_name}.smi in format: SMILES\\t ID\\t PREDICTION\\t DEVIATION\\n.")
    else:
        print(f"Gathered predictions were saved in {outdir}/{input_name}.smi in format: SMILES\\t ID\\t PREDICTION\\n.")

def write_results_sdf(g, input_name, outdir, label_id, original_sdf='', original_id_label='', original_prop_label='', prediction=False):
    """
    Writes the predictions made during optimization into an SDF file.
    @param g: igraph to be evaluated.
    @param input_name: Name of the input file to adapt the name of the written output.
    @param outdir: Name of the output directory.
    @param original_sdf: If provided, reads the molecules and writes its properties (label from argument original_prop_label) into the resulting SDF.
    @param original_id_label: Label to identify the original molecules.
    @param original_prop_label: Label containing the actual properties which should be stored in the resulting SDF.
    @param prediction: Decide whether there are vertex deviations annotated or not.
    """

    original_dict = None
    if original_sdf != '':
        original_dict = sdf_to_dict(original_sdf, original_id_label, original_prop_label)

    prop = f"{g['property']}_predicted"
    opt_prop_label = f"{prop}_optimized"

    w = Chem.SDWriter(f"{path.join(outdir, input_name)}.sdf")
    for vertex in g.vs.select(color='red'):
        mol = Chem.MolFromSmiles(vertex["name"])
        mol.SetProp(label_id, vertex["chem_id"])

        if original_dict is not None:
            or_entry = original_dict.get(vertex["chem_id"])
            if or_entry is None:
                print(f"Warning: Could not find molecule id {vertex['chem_id']} in {original_sdf}", file=sys.stderr)
                or_entry = ""
            mol.SetProp(original_prop_label, or_entry)

        mol.SetProp(prop, str(vertex["property"]))

        if prediction:
            mol.SetProp("deviation", str(vertex["variable"]))
            opt_prop = vertex["property"] + vertex["variable"]
            mol.SetProp(opt_prop_label, str(opt_prop))

        w.write(mol)
    print(f"Gathered predictions were saved in {path.join(outdir, input_name)}.sdf in SD format.")


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description="Optimizes an MMS network in GraphML format.",
                                     usage="python3 optimize_network.py FILE [-t, -s, -o, -l]")
    parser.add_argument('input_file', metavar='FILE.GraphML', type=str, help='Network to be optimized in GraphML format.')
    parser.add_argument('-t', '--transformations', dest='trafo_weight', default=0, type=float,
                        help="Decimal number, which indicates how strongly the number of transformations at an edge influence the sum of squares")
    parser.add_argument('-s', '--standarddeviation', dest='stddev', default=0, type=float,
                        help="Decimal number, which indicates how strongly the standard deviation of the transformation at an edge influence the sum of squares")
    parser.add_argument('-o', '--outdir', dest='outdir', default='optimized_results', type=str,
                        help='Output directory to store optimized network.')
    parser.add_argument('-l', '--logarithmic', dest='loga', default=False, action='store_true',
                        help='If this action is chosen, the number of transformations is considered logarithmically.')
    parser.add_argument('-f', '--format', metavar="FORMAT", choices=['sdf', 'smi'], default='smi', type=str,
                        help='Type of result file containing the predictions. Default is SMILES.')
    parser.add_argument('-i', '--sdf_id_label', default='_Name', type=str,
                        help='Name of the field to save the molecule id. Default is _Name (which is the name-header of a SDF entry).')
    parser.add_argument('-r', '--sdf_original', default='', type=str,
                        help='Original SD-file containing experimental data for the used test-set.')
    parser.add_argument('-x', '--sdf_original_property_label', default='', type=str,
                        help='Property label which is used of the original SD-file.')
    parser.add_argument('-y', '--sdf_original_id_label', default='', type=str,
                        help='ID label to identify compounds of the original SD-file.')
    parser.add_argument('-p', '--prediction', dest="prediction", action="store_true", default=False,
                        help="Optimize existing predictions: Create new fields deviation and '<property_name>_optimized'")
    args = parser.parse_args()

    if args.format == 'sdf':
        z = 1 if args.sdf_original != '' else 0
        z += 1 if args.sdf_original_property_label != '' else 0
        z += 1 if args.sdf_original_id_label != '' else 0
        if z > 0 and z < 3:
            print("Warning: If you need any experimental data from the original SDF you need to use all three arguments --sdf_original, --sdf_original_property_label, --sdf_original_id_label", file=sys.stderr)
            sys.exit(1)
    else:
        if args.sdf_original != '' or args.sdf_original_property_label != '' or args.sdf_original_id_label != '':
            print("Warning: You provided used the parameters sdf_original, sdf_original_property_label or sdf_original_id_label. These are not used when writing a smi file. Use SDF format if you need these information within the result file.", file=sys.stderr)


    # Set variables
    input_path = args.input_file
    input_name = path.splitext(path.basename(str(path.splitext(input_path)[0])))[0]
    pred = args.prediction
    try:
        # Construct igraph from given file
        g = ig.Graph.Read_GraphML(input_path)
    except FileNotFoundError:
        print(f"Your file was not found. No such file or directory: {input_path}")
        print("Please check your input and try again.")
        sys.exit(1)
    except:
        print("An error occurred. Maybe your file has the wrong format for this option.")
        print("Please upload a file in GraphML format, which was generated by a program from NetworkCreation.")
        sys.exit(1)

    # The actual algorithm to be performed
    trafo_weight = float(args.trafo_weight)
    stddev = float(args.stddev)
    optimize_network(g, trafo_weight, stddev, input_name, args.outdir, args.loga, pred)
    print(
        f"Note that the parameter following t is the weighting of edges corresponding to the number of transformations it belongs to ({trafo_weight})"
        f" and the parameter following s ist the influence of the standard deviation"
        f" ({stddev}) in the sum of squares.")

    if args.format == 'sdf':
        write_results_sdf(g, input_name, args.outdir, args.sdf_id_label, args.sdf_original, args.sdf_original_id_label, args.sdf_original_property_label, pred)
    else:
        write_results_smi(g, input_name, args.outdir, pred)
