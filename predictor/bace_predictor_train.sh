# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/bin/bash

training_sdf="data/bace_80_10_10/bace_train_tvt.sdf"
validation_sdf="data/bace_80_10_10/bace_validation_tvt.sdf"
test_sdf="data/bace_80_10_10/bace_test_tvt.sdf"
prop_label="pic50"
model_dir="predictor/models/bace_mpnn"
prediction_log="${model_dir}/prediction.log"
model="mpnn"

conda_env="predictor/nbs_dc_predictor_env.yml"
conda_env_name="nbs_dc_predictor_env"

# Activate conda environment
echo "Activating conda environment."
eval "$(conda shell.bash hook)"
env_check=$(conda env list | grep ${conda_env_name} | wc -l)
if (($env_check == 0)); then
  conda env create -f ${conda_env} -q
fi
conda activate ${conda_env_name}

# Build model
mkdir ${model_dir}
python3 predictor/build_model.py ${training_sdf} ${validation_sdf} ${test_sdf} ${prop_label} -m ${model} -o ${model_dir} | tee ${prediction_log}
conda deactivate
