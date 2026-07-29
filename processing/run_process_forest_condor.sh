#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo "Usage: $0 REPOSITORY EXECUTABLE INPUT_LIST OUTPUT_FILE MC_TYPE IS_PB_GOING PT_HAT TRIGGER_ID RECO_JET_SELECTION" >&2
    exit 2
fi

repository=$1
executable=$2
input_list=$3
output_file=$4
mc_type=$5
is_pb_going=$6
pt_hat=$7
trigger_id=$8
reco_jet_selection=$9

cd "$repository"
mkdir -p "$(dirname "$output_file")"

echo "Host:                  $(hostname)"
echo "Input list:            $input_list"
echo "Output file:           $output_file"
echo "MC type:               $mc_type"
echo "Pb-going:              $is_pb_going"
echo "pT-hat sample:         $pt_hat"
echo "Trigger ID:            $trigger_id"
echo "Reco-jet selection:    $reco_jet_selection"

exec "$executable" \
    "$input_list" \
    "$output_file" \
    "$mc_type" \
    "$is_pb_going" \
    "$pt_hat" \
    "$trigger_id" \
    "$reco_jet_selection"
