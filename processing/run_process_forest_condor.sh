#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 13 ]]; then
    echo "Usage: $0 REPOSITORY EXECUTABLE INPUT_LIST OUTPUT_FILE MC_TYPE IS_PB_GOING PT_HAT TRIGGER_ID RECO_JET_SELECTION VERTEX_FILTER_SELECTION COLLISION_SYSTEM JER_SMEARING_MODE USE_ETA_DEPENDENT_JER" >&2
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
vertex_filter_selection=${10}
collision_system=${11}
jer_smearing_mode=${12}
# Submitter defaults to 0: keep base standard and eta-dependent Extra families.
# Value 1 explicitly makes both families eta-dependent.
use_eta_dependent_jer=${13}

cd "$repository"
mkdir -p "$(dirname "$output_file")"

job_output=$(mktemp "${TMPDIR:-/tmp}/jetAnalysis-processForest.XXXXXX")
trap 'rm -f "$job_output"' EXIT

status=0
if "$executable" \
        "$input_list" \
        "$output_file" \
        "$mc_type" \
        "$is_pb_going" \
        "$pt_hat" \
        "$trigger_id" \
        "$reco_jet_selection" \
        "$vertex_filter_selection" \
        "$collision_system" \
        "$jer_smearing_mode" \
        "$use_eta_dependent_jer" >"$job_output" 2>&1; then
    :
else
    status=$?
fi

# processForestSimple can skip a transiently unreadable file and continue. A
# successful-looking multi-file output must not hide that data loss: retain a
# bounded diagnostic and let retry_failed_io_jobs.py resubmit the queue row.
read_error=0
if grep -Eq \
        'SysError in <TFile::ReadBuffer>:.*Input/output error|Error in <TBranch::GetBasket>:.*badread=|Failed to read chain entry .* from tree .* file ' \
        "$job_output"; then
    read_error=1
    status=74
fi

if [[ $status -ne 0 ]]; then
    {
        echo "processForestSimple failed (status $status)"
        echo "Host: $(hostname)"
        echo "Input list: $input_list"
        echo "Output file: $output_file"
        if [[ $read_error -eq 1 ]]; then
            echo "Detected ROOT read errors:"
            grep -E \
                'SysError in <TFile::ReadBuffer>:.*Input/output error|Error in <TBranch::GetBasket>:.*badread=|Failed to read chain entry .* from tree .* file ' \
                "$job_output" | tail -n 20
        fi
        echo "Last 64 KiB of combined application output:"
        tail -c 65536 "$job_output"
    } >&2
fi

exit "$status"
