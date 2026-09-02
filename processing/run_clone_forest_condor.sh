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

echo "Host:             $(hostname)"
echo "Input file:       $input_file"
echo "Output directory: $output_directory"
echo "Output file:      $output_file"

root -l -b -q "$macro_call"

# cloneForest reports failures from a void function, so ROOT can still return
# zero.  Make missing or empty output visible to Condor as a failed job.
if [[ ! -s "$output_file" ]]; then
    echo "Error: cloneForest did not create a nonempty output: $output_file" >&2
    exit 1
fi
