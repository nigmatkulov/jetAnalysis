#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 MACRO INPUT_FILE OUTPUT_DIRECTORY" >&2
    exit 2
fi

macro=$1
input_file=$2
output_directory=$3

mkdir -p "$output_directory"

# ROOT receives the macro arguments as a C++ expression.  The submitter rejects
# quotes and backslashes in paths, so interpolation here cannot change the
# expression structure.
macro_call="${macro}(\"${input_file}\",\"${output_directory}\")"
output_file="${output_directory%/}/${input_file##*/}"
job_output=$(mktemp "${TMPDIR:-/tmp}/jetAnalysis-cloneForest.XXXXXX")
trap 'rm -f "$job_output"' EXIT

status=0
if root -l -b -q "$macro_call" >"$job_output" 2>&1; then
    :
else
    status=$?
fi

# cloneForest reports failures from a void function, so ROOT can still return
# zero.  Make missing or empty output visible to Condor as a failed job.
if [[ ! -s "$output_file" ]]; then
    status=1
fi

if [[ $status -ne 0 ]]; then
    {
        echo "cloneForest failed (status $status)"
        echo "Host: $(hostname)"
        echo "Input file: $input_file"
        echo "Output file: $output_file"
        echo "Last 64 KiB of combined ROOT output:"
        tail -c 65536 "$job_output"
    } >&2
fi

exit "$status"
