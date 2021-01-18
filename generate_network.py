# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/usr/bin/env python3

import sys
import re
import igraph as ig
from statistics import mean, stdev
import argparse
from collections import defaultdict

def read_indexing_file(file_name, property_dict, property, predictions={}):
    """
    Processing of the indexing file provided from the output of mmpdb.
    :param file_name: Path to the indexing file
    :param property_dict: Dictionary of molecules with known properties as molecule IDs as keys and corresponding
                          properties as values
    :param property: The property type of interest (e.g. pic50)
    :param predictions: Dictionary of molecules with existing predictions for their properties represented by
                        molecule IDs as keys and corresponding property predictions as values
    :return: nodes, edges, transformations and ids as dictionaries (see comments in code)
    """
    # Node key: SMILES string, value: properties for all occurrences of the same SMILES to calculate the mean for this
    # SMILES expression.
    nodes = defaultdict(list)
    # Edge key: node pair (a,b), value: transformation
    edges = defaultdict(list)
    # Transformation key: SMIRKS expression, value: List of all property changes associated with this transformation
    transformations = {}
    # ID key: SMILES, value: ID associated with this smiles.
    ids = {}
    # Counter for red vertices with unknown and unpredicted property
    nof_unknown_props = 0
    print("Constructing igraph, this may take some time.")
    for num, line in enumerate(open(file_name, 'r')):
        line_fields = re.split('\s|,', line)
        # Indexing file 'csv' format is a tab-separated table with the columns:
        #   SMILES1  SMILES2  id1  id2  V1>>V2  C
        left_smi = line_fields[0]
        right_smi = line_fields[1]
        left_id = line_fields[2]
        right_id = line_fields[3]
        if left_id in predictions.keys():
            left_property = predictions[left_id][property]
        else:
            left_property = property_dict[left_id][property]
        if right_id in predictions.keys():
            right_property = predictions[right_id][property]
        else:
            right_property = property_dict[right_id][property]
        transformation = line_fields[4]
        ids[left_smi] = left_id
        ids[right_smi] = right_id

        # Collecting node set for igraph
        try:
            nodes[left_smi].append(float(left_property))
        except ValueError:
            nof_unknown_props += 1
            left_property = '?'
        try:
            nodes[right_smi].append(float(right_property))
        except ValueError:
            nof_unknown_props += 1
            right_property = '?'

        # Collecting edge set for igraph.
        try:
            edges[left_smi, right_smi].append(transformation)
        except ValueError:
            print(f"Invalid transformation: {transformation} in line {num}")

        # Collecting information about transformations for mean calculation:
        if not transformation in transformations:
            transformations[transformation] = []
        if left_property != '?' and right_property != '?':
            trafo_property = float(right_property) - float(
                left_property)
            transformations[transformation].append(trafo_property)
    if nof_unknown_props > 0 :
        if len(predictions) > 0:
            print("WARNING: Prediction optimization mode was selected, but molecules with invalid or unexisting predictions were found. Their property was treated as unknwon.")
        print(f"Number of vertices with unknown property: {nof_unknown_props}")
    # Return all collected information
    return nodes, edges, transformations, ids

def calc_means(transformations):
    """
    For each transformation, calculate the mean of its associated property changes.
    @param transformations: Dictionary containing a list of all property changes occurring with a transformation
    @return: The mean of the corresponding transformations as dictionary,
             the number of values for each transformation and the standard deviation of its property change occurrences.
    """
    num_vals = {}
    std_dev = {}
    for trafo in transformations:
        num_vals[trafo] = len(transformations[trafo])
        if len(transformations[trafo]) > 1:
            std_dev[trafo] = stdev(transformations[trafo])
        else:
            std_dev[trafo] = 0
        if transformations[trafo] != []:
            transformations[trafo] = mean(transformations[trafo])
        else:
            transformations[trafo] = 0.0  # We assume there is no change for unknown transformations
    return transformations, num_vals, std_dev


def construct_igraph(file_name, property_dict, property, predictions={}):
    """
    Constructs an igraph object from a given mmpdb indexing csv output.
    Nodes and edges save mean values as their properties.
    Edges are annotated by the mean property change between matched molecular pairs that are connected by the same transformation.
    @param file_name: The input file taken from mmpdb indexing output
    @param property_dict: Dictionary of molecules with known properties as molecule IDs as keys and corresponding
                          properties as values
    @param property: The property type of interest (e.g. pic50)
    @param predictions: Dictionary of molecules with existing predictions for their properties represented by
                        molecule IDs as keys and corresponding property predictions as values
    @return: An igraph object with following anntoation:
            - at the nodes: name (equals the SMILES string), property (the mean of all ChEMBL IDs with the same SMILES)
            - at the edges: transformation, transform_mean (the mean of all property changes between MMP molecules connected by this transformation),
                transform_num (number of occurrences of the transformation), std_dev (standard deviation belonging to transform_mean)
    """
    try:
        nodes, edges, transformations, ids = read_indexing_file(file_name, property_dict, property, predictions)
    except KeyError:
        print("ERROR: Your indexing file contains an ID which is not present in your input smiles file(s).")
        print("Please make sure that your indexing file was created from your input smiles file.")
        sys.exit(1)
    print("Done reading indexing file, constructing graph.")
    # Construct igraph from dictionary information:
    g = ig.Graph(directed=True)
    pred_key_set = set(predictions.keys())
    for node in nodes:
        # Color vertices by their property: green = known property, red = unknown property
        if nodes[node] != []:
            if ids[node] in pred_key_set:
                g.add_vertex(name=node, property=mean(nodes[node]), color='red', chem_id=ids[node])
            else:
                g.add_vertex(name=node, property=mean(nodes[node]), color='green', chem_id=ids[node])
        else:
            g.add_vertex(name=node, property=0.0, color='red', chem_id=ids[node])

    transformations, num_vals, std_dev = calc_means(transformations)

    # Construct the edges, care about each mmp having only two antiparallel edges connecting them.
    black_list = set()
    for num, edge in enumerate(edges):
        if num % 1000 == 0 or num == len(edges) - 1:
            text = f"\rConstructing edge: {1 + num}/{len(edges)} ({float('{:.2f}'.format(((1 + num) / len(edges)) * 100))})%"
            sys.stdout.write(text)
            sys.stdout.flush()
        if edge not in black_list:
            trafo_list = []
            for trafo in edges[edge]:
                trafo_list.append((trafo, transformations[trafo]))
            chosen_trafo = max(trafo_list, key=lambda item: item[1])
            g.add_edge(edge[0], edge[1], transformation=chosen_trafo[0], transform_mean=chosen_trafo[1],
                       std_dev=std_dev[chosen_trafo[0]], transform_num=num_vals[chosen_trafo[0]])
            g.add_edge(edge[1], edge[0], transformation=chosen_trafo[0].replace(">>", "<<"),
                       transform_mean=-chosen_trafo[1], std_dev=std_dev[chosen_trafo[0]],
                       transform_num=num_vals[chosen_trafo[0]])
            black_list.add((edge[1], edge[0]))

    for edge in g.es:
        edge["label"] = str(edge['transform_mean'])

    for vertex in g.vs:
        if vertex["color"] == 'green':
            vertex["label"] = f"{vertex['name']} \n  {property}: {vertex['property']}"
        else:
            vertex["label"] = f"{vertex['name']} \n  {property}: ?"
    print()
    print("Done constructing graph.")
    g["property"] = property
    return g


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Create an MMS network in .GraphML format',
                                     usage="python3 generate_network.py INDEX.csv MOLECULES.sdf PROPERTY [-o]")
    parser.add_argument('index', metavar="INDEX.csv", type=str,
                        help='Indexed output from mmpdb containing the MMS relations.')
    parser.add_argument('smi_input', metavar="INPUT.smi", type=str,
                        help='Input file in smiles format containing the molecules. Row format: {SMILES\\t ID\\t PROPERTY}.')
    parser.add_argument("property", metavar="PROPERTY", type=str,
                        help='Property of interest, annotated in the input SD file.')
    parser.add_argument('-o', '--outname', metavar='OUTNAME', default='graph.GraphML',
                        help='Specify name of output graph file.')
    parser.add_argument('-p', '--prediction', metavar="PREDICTION.smi", type=str,
                        help="Predictions in smiles format of INPUT.smi")
    args = parser.parse_args()

    smi_file = args.smi_input # Input molecules
    idx_file = args.index # mmpdb result (indexed file)
    property = args.property # Property of interest
    outname = args.outname # Name of the output result file
    molecules = {} # All molecules associated with their properties

    # Extraction of molecules, their IDs and properties
    for line in open(smi_file):
        infos = line.strip().split('\t')  # idx: 0 = smiles, 1 = id, 2 = property
        molecules[infos[1]] = {'smiles': infos[0], property: infos[2] if infos[2] != "nan" else ""}

    # Construction of igraph with or without already existing predictions from another predictor model
    if args.prediction:
        pred_dict = {}
        for line in open(args.prediction):
            infos_pred = line.strip().split('\t')  # idx: 0 = smiles, 1 = id, 2 = property
            pred_dict[infos_pred[1]] = {'smiles': infos_pred[0], property: infos_pred[2] if infos_pred[2] != "nan" else ""}
        g = construct_igraph(idx_file, molecules, property, pred_dict)
    else:
        g = construct_igraph(idx_file, molecules, property)
    g.write_graphml(f"{outname}")
    print(f"Your network was saved in GraphML format as: {outname}")
