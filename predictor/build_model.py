# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey
 
# Parts of the DeepChem library (https://github.com/deepchem/) are included.
# Thus, this script is additionally licensed unter MIT License.

#!/usr/bin/env python3

import numpy as np
import deepchem
from deepchem.molnet.preset_hyper_parameters import hps
from deepchem.utils.save import save_dataset_to_disk
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
import argparse


def build_model(training_sdf, validation_sdf, test_sdf, featurizer_func, model, model_dir, features):
    loader = deepchem.data.SDFLoader(tasks=[args.task], clean_mols=True, featurizer=featurizer_func)
    print(loader)
    train_dataset = loader.featurize(training_sdf)
    print(f"train has {len(train_dataset)} compounds")
    valid_dataset = loader.featurize(validation_sdf)
    test_dataset = loader.featurize(test_sdf)

    transformers = [
        deepchem.trans.NormalizationTransformer(
            transform_y=True, dataset=train_dataset, move_mean=True)
    ]

    for transformer in transformers:
        train_dataset = transformer.transform(train_dataset)
        valid_dataset = transformer.transform(valid_dataset)
        test_dataset = transformer.transform(test_dataset)

    tasks = [args.task]
    n_features = features
    metric = [deepchem.metrics.Metric(
        deepchem.metrics.rms_score, np.mean,
        mode='regression'),
        deepchem.metrics.Metric(
            deepchem.metrics.pearson_r2_score,
            np.mean,
            mode='regression')]
    train_scores, valid_scores, test_scores = benchmark_regression(train_dataset,
                                                                   valid_dataset,
                                                                   test_dataset,
                                                                   tasks,
                                                                   transformers,
                                                                   n_features,
                                                                   metric,
                                                                   model,
                                                                   test=True,
                                                                   model_dir=model_dir)
    print("Train score:", train_scores)
    print("Valid score:", valid_scores)
    print("Test score:", test_scores)
    print("Model was saved to:", model_dir)


# This function was copied and modified after the existing example of deepchem
def benchmark_regression(train_dataset,
                         valid_dataset,
                         test_dataset,
                         tasks,
                         transformers,
                         n_features,
                         metric,
                         model,
                         test=False,
                         hyper_parameters=None,
                         seed=123,
                         model_dir='model_dir'):
    """
    Calculate performance of different models on the specific dataset & tasks

    Parameters
    ----------
    train_dataset: dataset struct
        dataset used for model training and evaluation
    valid_dataset: dataset struct
        dataset only used for model evaluation (and hyperparameter tuning)
    test_dataset: dataset struct
        dataset only used for model evaluation
    tasks: list of string
        list of targets(tasks, datasets)
    transformers: dc.trans.Transformer struct
        transformer used for model evaluation
    n_features: integer
        number of features, or length of binary fingerprints
    metric: list of dc.metrics.Metric objects
        metrics used for evaluation
    model: string, optional
        choice of model
        'tf_regression', 'tf_regression_ft', 'rf_regression', 'graphconvreg',
        'dtnn', 'dag_regression', 'xgb_regression', 'weave_regression',
        'textcnn_regression', 'krr', 'ani', 'krr_ft', 'mpnn'
    test: boolean, optional
        whether to calculate test_set performance
    hyper_parameters: dict, optional (default=None)
        hyper parameters for designated model, None = use preset values
    model_dir: directory to save the model


    Returns
    -------
    train_scores : dict
      predicting results(R2) on training set
    valid_scores : dict
      predicting results(R2) on valid set
    test_scores : dict
      predicting results(R2) on test set

    """
    train_scores = {}
    valid_scores = {}
    test_scores = {}
    assert model in ['graphconvreg', 'mpnn']
    if hyper_parameters is None:
        hyper_parameters = hps[model]
    model_name = model

    if model_name == 'graphconvreg':
        batch_size = hyper_parameters['batch_size']
        nb_epoch = hyper_parameters['nb_epoch']
        learning_rate = hyper_parameters['learning_rate']
        n_filters = hyper_parameters['n_filters']
        n_fully_connected_nodes = hyper_parameters['n_fully_connected_nodes']

        model = deepchem.models.GraphConvModel(
            len(tasks),
            graph_conv_layers=[n_filters] * 2,
            dense_layer_size=n_fully_connected_nodes,
            batch_size=batch_size,
            learning_rate=learning_rate,
            random_seed=seed,
            mode='regression',
            model_dir=model_dir)

    elif model_name == 'mpnn':
        batch_size = hyper_parameters['batch_size']
        nb_epoch = hyper_parameters['nb_epoch']
        learning_rate = hyper_parameters['learning_rate']
        T = hyper_parameters['T']
        M = hyper_parameters['M']

        model = deepchem.models.MPNNModel(
            len(tasks),
            n_atom_feat=n_features[0],
            n_pair_feat=n_features[1],
            n_hidden=n_features[0],
            batch_size=batch_size,
            learning_rate=learning_rate,
            use_queue=False,
            mode="regression",
            model_dir=model_dir
        )

    print('-----------------------------')
    print('Start fitting: %s' % model_name)
    if nb_epoch is None:
        model.fit(train_dataset)
    else:
        model.fit(train_dataset, nb_epoch=nb_epoch)
    model.restore()
    save_dataset_to_disk(model_dir, train_dataset, valid_dataset, test_dataset, transformers)
    print('Evaluating')
    train_scores[model_name] = model.evaluate(train_dataset, metric, transformers)
    valid_scores[model_name] = model.evaluate(valid_dataset, metric, transformers)
    if test:
        test_scores[model_name] = model.evaluate(test_dataset, metric, transformers)
    return train_scores, valid_scores, test_scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict for a data split into training, validation and test data.",
                                     usage="predictor.py")
    parser.add_argument("training_SDF", metavar="TRAINING_SDF", type=str, help="Training set in SD format.")
    parser.add_argument("validation_SDF", metavar="VALIDATION_SDF", type=str, help="Validation set in SD format.")
    parser.add_argument("test_SDF", metavar="TEST_SDF", type=str, help="Test set in SD format.")
    parser.add_argument("task", metavar="TASK", type=str,
                        help="Task for the model to perform. Depends on the input data set.",
                        choices=["pic50", "logD", "expt", "sol"])
    parser.add_argument("-o", "--outdir", metavar="OUTDIR", type=str, default="model_dir",
                        help="Directory to save the model")
    parser.add_argument("-m", "--model", metavar="MODEL", type=str, choices=["graphconvreg", 'mpnn'],
                        help="Model type to be constructed.")
    args = parser.parse_args()

    model_type = args.model
    model_dir = args.outdir
    train = args.training_SDF
    valid = args.validation_SDF
    test = args.test_SDF

    if model_type == "graphconvreg":
        featurizer_func = deepchem.feat.ConvMolFeaturizer()
        print("Starting Graph Convolution for Dataset")
        build_model(train, valid, test, featurizer_func, model_type, model_dir, 75)
    elif model_type == "mpnn":
        featurizer_func = deepchem.feat.WeaveFeaturizer()
        print("Starting MPNN for Dataset")
        build_model(train, valid, test, featurizer_func, model_type, model_dir, [75, 14])
    else:
        print("Unsupported model choice. Please check your input arguments.")
