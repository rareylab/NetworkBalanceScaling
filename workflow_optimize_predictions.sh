#!/bin/bash

# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

source ~/.bashrc

training_sdf="data/bace_80_10_10/bace_train_tvt.sdf"
validation_sdf="data/bace_80_10_10/bace_validation_tvt.sdf"
query_sdf="data/bace_80_10_10/bace_test_tvt.sdf"
original="data/bace.sdf"

id_label="id"
prop_label="pic50"

temp_dir="workflow_optimize_predictions_tmp"
results_dir="bace_optimized_pred"
result_name="bace_result"

predictor="./predictor/bace_predictor_predict.sh"



training_sdf=`realpath $training_sdf`
query_sdf=`realpath $query_sdf`
predictor=`realpath $predictor`

merged_sdf="${temp_dir}/merged.sdf"
merged_smi="${temp_dir}/merged.smi"
graph_train_sdf="${temp_dir}/graph_train.sdf"
graph_train_smi="${temp_dir}/graph_train.smi"
merged_fragments="${temp_dir}/merged.fragments.gz"
merged_index="${temp_dir}/merged_index.csv"
mmpdb_idx_csv="${temp_dir}/merged_idx.csv"
graph="${temp_dir}/${result_name}.GraphML"

prediction_res="${results_dir}/prediction.smi"
prediction_log="${results_dir}/prediction.log"

evaluation_log="${results_dir}/evaluation.log"

conda_env="network_balance_scaling_env.yml"
conda_env_name="network_balance_scaling_env"

# 1. call predictor on query.sdf and save result in prediction_res_tmp
mkdir ${temp_dir}
mkdir ${results_dir}
echo "Creating prediction with deepchem model. This may take a while..."
# Predictor needs to create an output in smiles format, saved in prediction_res (prediction.smi)
# Note: The result of the prediction needs to follow the format: SMILES\tID\tPREDICTION\n
cmd="${predictor} ${query_sdf} ${prediction_res} | tee ${prediction_log}"
echo "run: ${cmd}"
eval ${cmd}

## 2. activate conda environment
echo "Activating conda environment."
eval "$(conda shell.bash hook)"
env_check=$(conda env list | grep ${conda_env_name} | wc -l)
if (($env_check == 0)); then
  conda env create -f ${conda_env} -q
fi
conda activate ${conda_env_name}

### 3. create merged sdf and smiles and call mmpdb
cat ${training_sdf} ${validation_sdf} ${query_sdf} > ${merged_sdf}
cat ${training_sdf} ${validation_sdf} > ${graph_train_sdf}
echo "Converting files with utils scripts."
python3 utils/convert_SDF2SMI.py ${merged_sdf} ${prop_label} --output ${merged_smi} --id ${id_label}
python3 utils/convert_SDF2SMI.py ${graph_train_sdf} ${prop_label} --output ${graph_train_smi} --id ${id_label}
python3 utils/mmpdb/mmpdb fragment ${merged_smi} -o ${merged_fragments}
python3 utils/mmpdb/mmpdb index ${merged_fragments} -o ${merged_index}

# 3. generate network
cmd="python3 generate_network.py ${merged_index} ${graph_train_smi} ${prop_label} -o ${graph} --prediction ${prediction_res}"
echo "run: ${cmd}"
eval ${cmd}

# 4. optimize network and existing predictions with standard parameters
cmd="python3 optimize_network.py ${graph} --outdir ${results_dir} --prediction --format sdf --sdf_id_label ${id_label} --sdf_original ${original} --sdf_original_id_label ${id_label} --sdf_original_property_label ${prop_label}"
echo "run: ${cmd}"
eval ${cmd}


# 5. evaluate optimized results
if [ "$original" != "" ]; then
    label_predicted="${prop_label}_predicted"
    label_optimized="${label_predicted}_optimized"
    cmd="python3 utils/evaluate_results.py -o ${original} -p ${prediction_res} -q ${results_dir}/${result_name}.sdf -i ${id_label} -l ${evaluation_log} --label_predicted ${label_predicted} --label_optimized ${label_optimized} ${prop_label}"
    echo "run ${cmd}"
    eval ${cmd}
fi

##### 6. delete temp dir & deactivate environment
result_file="${results_dir}/${result_name}_optimized.GraphML"
if [ -f ${result_file} ] && [ -s ${result_file} ]; then
  # rm -r ${temp_dir}
    # echo "Deleted temporary files. Results can be found in ${results_dir}."
  echo "Ok."
else
  echo "An error occurred. Check temporary results in ${temp_dir}."
fi
conda deactivate
