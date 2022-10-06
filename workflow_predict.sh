#!/bin/bash

# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey


training_sdf="data/bace_train.sdf"
query_sdf="data/bace_test.sdf"
original="data/bace.sdf"

id_label="id"
id_label_prediction="_Name"
prop_label="pic50"

temp_dir="workflow_predict_tmp"

merged_sdf="${temp_dir}/merged.sdf"
merged_smi="${temp_dir}/merged.smi"
merged_fragments="${temp_dir}/merged.fragments.gz"
merged_index="${temp_dir}/merged_index.csv"
mmpdb_idx_csv="${temp_dir}/merged_idx.csv"
result_name="bace_split"

graph="${temp_dir}/${result_name}.GraphML"
results_dir="bace_optimized"
evaluation_log="${results_dir}/evaluation.log"

conda_env="network_balance_scaling_env.yml"
conda_env_name="network_balance_scaling_env"


# 1. activate conda environment
echo "Activating conda environment."
eval "$(conda shell.bash hook)"
env_check=$(conda env list | grep ${conda_env_name} | wc -l)
if (($env_check == 0)); then
  conda env create -f ${conda_env} -q
fi
conda activate ${conda_env_name}

# 2. create merged sdf and smiles and call mmpdb
echo "Creating mmpdb index."
mkdir ${temp_dir}
cat ${training_sdf} ${query_sdf} > ${merged_sdf}
cmd="python3 utils/convert_SDF2SMI.py ${merged_sdf} ${prop_label} --output ${merged_smi} --id ${id_label}"
echo "run: ${cmd}"
eval ${cmd}

cmd="python3 utils/mmpdb/mmpdb fragment ${merged_smi} -o ${merged_fragments}"
echo "run: ${cmd}"
eval ${cmd}

cmd="python3 utils/mmpdb/mmpdb index ${merged_fragments} -o ${merged_index}"
echo "run: ${cmd}"
eval ${cmd}

# 3. generate network
echo "Generating network."

cmd="python3 generate_network.py ${merged_index} ${merged_smi} ${prop_label} -o ${graph}"
echo "run: ${cmd}"
eval ${cmd}

# 4. optimize network with standard parameters
echo "Optimizing network."
cmd="python3 optimize_network.py ${graph} --outdir ${results_dir} --format sdf"
echo "run: ${cmd}"
eval ${cmd}

# 5. evaluate optimized results

label_predicted="${prop_label}_predicted"
label_optimized="${label_predicted}_optimized"
cmd="python3 utils/evaluate_results.py -o ${original} -p ${results_dir}/${result_name}.sdf --label_predicted ${label_predicted} --label_optimized ${label_optimized} -i ${id_label} --id_prediction ${id_label_prediction} -l ${evaluation_log} ${prop_label}"
echo "run: ${cmd}"
eval ${cmd}

# 6. delete temp dir & deactivate environment
result_file="${results_dir}/${result_name}_optimized.GraphML"
if [ -f ${result_file} ] && [ -s ${result_file} ]; then
  rm -r ${temp_dir}
  echo "Deleted temporary files. Results can be found in ${results_dir}."
else
  echo "An error occurred. Check temporary results in ${temp_dir}."
fi

conda deactivate
