# This file is part of Network Balance Scaling, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2021
# University of Hamburg, ZBH - Center for Bioinformatics
# Sophia Hönig, Wolf-Guido Bolick, Emanuel Ehmki, Matthias Rarey

#!/bin/bash

dirname=`dirname $0`
script_dir=`cd $dirname && pwd`
templ_file="${script_dir}/template_header_opt_pred.txt"

complete_workflow="${script_dir}/workflow_optimize_predictions.sh"


help="Syntax: ./create_optimize_pred_workflow.sh <train_sdf> <query_sdf> <id_label> <prop_label> <result_name> <result_dir> <predictor_script> [validation_sdf] [original_sdf]"

if [ "$#" -lt 6 ]; then
    echo "Illegal number of parameters" >> /dev/stderr
    echo "$help"
    exit -1
fi

tr_sdf="$1"
qy_sdf="$2"
id_lab="$3"
pr_lab="$4"
rs_nam="$5"
rs_dir="$6"
prd_sh="$7"
vl_sdf="${8:-/dev/null}"
or_sdf="${9:-}"


tmp_dir="tmp_${rs_dir}"
new_script_name="optimize_workflow_${rs_nam}.sh"

echo "train_sdf: $tr_sdf"
echo "query_sdf: $qy_sdf"
echo "id_label: $id_lab"
echo "prop_label: $pr_lab"
echo "result_name: $rs_nam"
echo "result_dir: $rs_dir"
echo "predictor_script: $prd_sh"
echo "validation_sdf: $vl_sdf"
echo "original_sdf: $or_sdf"


cat $templ_file | \
    sed "s%###_training_sdf_###%${tr_sdf}%" | \
    sed "s%###_validation_sdf_###%${vl_sdf}%" | \
    sed "s%###_query_sdf_###%${qy_sdf}%" | \
    sed "s%###_original_sdf_###%${or_sdf}%" | \
    sed "s%###_id_label_###%${id_lab}%" | \
    sed "s%###_prop_label_###%${pr_lab}%" | \
    sed "s%###_temp_dir_###%${tmp_dir}%" | \
    sed "s%###_results_dir_###%${rs_dir}%" | \
    sed "s%###_result_name_###%${rs_nam}%" | \
    sed "s%###_predictor_sh_###%${prd_sh}%" > $new_script_name

ns_lines=`wc -l $new_script_name | cut -d ' ' -f 1`


cmd="cat ${complete_workflow} | tail -n +${ns_lines} >> $new_script_name; chmod +x $new_script_name"

echo "cmd is $cmd"
eval $cmd
