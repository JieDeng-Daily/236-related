#!/usr/bin/bash

folder_with_pdbs=${1:-"./inputs/"}
output_dir=${2:-"./outputs/"}
chains_to_design=${3:-"A"}
design_only_positions=${4:-"5 413 174 177"}

if [ "$1" == "-h" ]; then
    echo "Usage: mpnn_design_specific.sh [folder_with_pdbs] [output_dir] [chains_to_design] [design_only_positions]"
    echo "Defaults:"
    echo "  folder_with_pdbs: ./inputs/"
    echo "  output_dir: ./outputs/"
    echo "  chains_to_design: A"
    echo "  design_only_positions: 5 413 174 177"
    exit 0
fi

if [ ! -d $folder_with_pdbs ]; then
    echo "Error: Input folder $folder_with_pdbs does not exist."
    exit 1
fi

if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"

echo "Starting ProteinMPNN pipeline..."
echo "Input folder: $folder_with_pdbs"
echo "Output directory: $output_dir"
echo "Chains to design: $chains_to_design"
echo "Design only positions: $design_only_positions"

python /data/dj/softwares/ProteinMPNN-main/helper_scripts/parse_multiple_chains.py \
    --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
if [ $? -ne 0 ]; then
    echo "Error occurred in parsing chains step."
    exit 1
fi

python /data/dj/softwares/ProteinMPNN-main/helper_scripts/assign_fixed_chains.py \
    --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
if [ $? -ne 0 ]; then
    echo "Error occurred in assigning fixed chains step."
    exit 1
fi

python /data/dj/softwares/ProteinMPNN-main/helper_scripts/make_fixed_positions_dict.py \
    --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions \
    --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed
if [ $? -ne 0 ]; then
    echo "Error occurred in creating fixed positions dictionary step."
    exit 1
fi

python /data/dj/softwares/ProteinMPNN-main/protein_mpnn_run.py \
    --jsonl_path $path_for_parsed_chains \
    --chain_id_jsonl $path_for_assigned_chains \
    --fixed_positions_jsonl $path_for_fixed_positions \
    --out_folder $output_dir \
    --num_seq_per_target 1000 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1
if [ $? -ne 0 ]; then
    echo "Error occurred in ProteinMPNN run step."
    exit 1
fi

echo "ProteinMPNN pipeline completed successfully."

