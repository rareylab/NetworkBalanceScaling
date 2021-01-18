# Usage
## Workflows
Please execute the following workflows depending on your preferences:
### workflow_predict.sh
- installs the required conda environment if necessary
- generates a network
- optimizes the network and predicts missing molecular properties
- output contains the predictions and the graph in separate files

### workflow_optimize_predictions.py
- installs the required conda environment if necessary
- predicts missing molecular properties using a specified machine learning script  
- generates a network containing the formerly made predictions
- optimizes the network including the annotated predictions 
- output contains the predictions with optimization suggestions and the graph in separate files

### Further Advices
- simply execute the shell scripts provided in this directory
- for the usage of Network Balance Scaling as a predictor, please try and modify workflow_predict.sh
- for the usage of Network Balance Scaling to optimize existing predictions, please try and modify workflow_optimize_predictions.sh
  - this workflow contains a precalculated model
  - new models can be created with our inhouse predictor based on deepchem using the script predidictors/bace_predictor_train.sh
  - to generate predictions with an existing model independent from the workflow, please try the script predictors/bace_predictor_predict.sh
- please note, that our workflows only include the BACE data set of moleculenet due to a huge RAM requirement for larger data sets of this benchmark
  - if you are interested in the other data sets of this benchmark, feel free to construct your own pipeline and try our helper scripts in the utils directory
  - you can find precalculated benchmark scenarios and the original CSV files from the moleculenet benchmark in the data directory
- **note** that all scripts should be startet from main directory in order of the path management
- to evaluate the results you achieved, please checkout our script **utils/evaluate_RMSD.py**


## Dependencies
- python is needed to run the scripts
- if you use an anaconda environment:
  - the workflows will automatically create and activate the required dependencies in extra environments
  - if the activation fails, please add the following line insinde the workflows above the activate command:
    - _source ~/anaconda3/etc/profile.d/conda.sh_
  - or add this line to your .bashrc
- if you use another python environment, please make sure to provide the below listed dependencies

### Tool Dependencies
- python-igraph=0.8.3=py37h340e831_2
- cvxpy=1.1.10=py37h89c1867_0
- cvxpy-base=1.1.10=py37hdc94413_0
- rdkit=2020.09.4=py37he53b9e1_0

## generate_network.py
A script for creating networks with or without existing predictions at red vertices.
```
usage: python3 generate_network.py INDEX.csv MOLECULES.sdf PROPERTY [-o]

Create an MMS network in .GraphML format

positional arguments:
  INDEX.csv             Indexed output from mmpdb containing the MMS
                        relations.
  INPUT.smi             Input file in smiles format containing the molecules.
                        Row format: {SMILES\t ID\t PROPERTY}.
  PROPERTY              Property of interest, annotated in the input SD file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTNAME, --outname OUTNAME
                        Specify name of output graph file.
  -p PREDICTION.smi, --prediction PREDICTION.smi
                        Predictions in smiles format of INPUT.smi

```
## optimize_network.py
A script for optimizing networks. It can be used to either predict unknwon properties or to optimize existing
predictions depending on the option ```-p```
```
usage: python3 optimize_network.py FILE [-t, -s, -o, -l]

Optimizes an MMS network in GraphML format.

positional arguments:
  FILE.GraphML          Network to be optimized in GraphML format.

optional arguments:
  -h, --help            show this help message and exit
  -t TRAFO_WEIGHT, --traformations TRAFO_WEIGHT
                        Decimal number, which indicates how strongly the
                        number of transformations at an edge influence the sum
                        of squares
  -s STDDEV, --standarddeviation STDDEV
                        Decimal number, which indicates how strongly the
                        standard deviation of the transformation at an edge
                        influence the sum of squares
  -o OUTDIR, --outdir OUTDIR
                        Output directory to store optimized network.
  -l, --logarithmic     If this action is chosen, the number of
                        transformations is considered logarithmically.
  -f FORMAT, --format FORMAT
                        Type of result file containing the predictions.
                        Default is SMILES.
  -p, --prediction      Optimize existing predictions
```

# Predictors

## Integrated Predictor
The integrated predictor is based on DeepChem (https://github.com/deepchem/).
Thus, the corresponding python scripts **predictor/build_model.py** and **predictor/predictor.py** additionally underly the MIT License. Our integrated predictor comprises:
- a trained graph convolution regression model for the BACE data set
- a trained message passing neural network model for the BACE data set
- the Python script **predictor/build_model.py** which can be used to create a model
- the Python script **predictor/predictor.py** which can be used to create predictions using an existing model
- the bash script **predictor/bace_predictor_train.sh** which trains a model based on the BACE data set
- the bash script **predictor/bace_predictor_predict.sh** which predicts properties for the BACE test set
- training, validation and test set for BACE can be found in our **data**
- for the other data sets, please use the scripts provided in **utils** to create your own data splits and convert files with respective formats

### Integrated Predictor Dependencies
- deepchem=2.3.0=py_2
- rdkit=2020.09.2=py37h713bca6_0
- tensorflow=1.14.0=h4531e10_0

## Integrating Your Own Predictor
To incorporate existing predictions into your network, please mind the following:
- the resulting predictions need to be written to a SMILES file
- the SMILES file needs to follow the format **SMILES**\t**ID**\t**PROPERTY**
- using **workflow_optimize_predictions.sh** you need to replace the call of our default predictor based on deepchem (_./predictor/bace_predictor_predict.sh_)
- please remember to adapt the training, validation and query set at the beginning of the bash script
- **note:** Please remember, that all molecules need to be included into the mmpdb index. Using the network generation with existing predictions independently from **workflow_optimize_predictions** 
  requires the molecules with known properties (training and validation sets) in a first SMILES file, your predictions in a second SMILES file and the index including all (training, validation and query) molecules.

# File Formats
- the workflows are constructed to run with SD files as input
- the SD files are converted during the workflow using the scripts provided at **utils/**
- for directly using SMILES files, please note that the required format consists of three tab-separated columns:
  - **SMILES**\t**ID**\t**PROPERTY**
- the SMILES format is also a requirement for the input files used in mmpdb's MMP analysis
