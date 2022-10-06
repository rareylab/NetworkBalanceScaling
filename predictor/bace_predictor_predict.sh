#!/bin/bash

# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey


conda_env="predictor/nbs_dc_predictor_env.yml"
conda_env_name="nbs_dc_predictor_env"
property="pic50"
model_dir="predictor/models/bace_mpnn"
model="mpnn"

default_query="data/bace_80_10_10/bace_test_tvt.sdf"
default_result="prediction.smi"

query_sdf=${1:-$default_query}
result=${2:-$default_result}

# Activate conda environment
echo "Activating conda environment."
eval "$(conda shell.bash hook)"
env_check=$(conda env list | grep ${conda_env_name} | wc -l)
if (($env_check == 0)); then
  conda env create -f ${conda_env} -q
fi
conda activate ${conda_env_name}

# Predict
cmd="python3 predictor/predictor.py ${query_sdf} ${property} -d ${model_dir} -o ${result} -m ${model} -i id"
echo "run: ${cmd}"
eval ${cmd}
conda deactivate
