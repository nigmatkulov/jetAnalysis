#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_file> <N> <sample_name>"
    echo "sample_name: 0 for DATA, 1 for PYTHIA"
    exit 1
fi

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Input file name
input_file="$1"
# Number of entries for each sublist
N="$2"

# Sample name (DATA or PYTHIA)
if [ "$3" -eq 0 ]; then
    sample_name="DATA"
    sample_prefix="pp5020_data"
else
    sample_name="PYTHIA"
    sample_prefix="pp5020_pythia"
fi

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file '$input_file' not found."
    exit 1
fi

# Count the total number of lines in the input file
total_lines=$(wc -l < "$input_file")

# Calculate the number of sublists needed
num_sublists=$((total_lines / N))
if [ $((total_lines % N)) -ne 0 ]; then
    ((num_sublists++))
fi

if [ ! -d "$PWD/input/pp5020/${formatted_date}" ]; then
    #echo "Directory '$PWD/input/pp5020/${formatted_date}' does not exist. Creating..."
    mkdir -p "$PWD/input/pp5020/${formatted_date}"
fi

# Create sublists
for ((i = 0; i < num_sublists; i++)); do
    start=$((i * N + 1))  # Calculate start line number for current sublist
    end=$((start + N - 1))  # Calculate end line number for current sublist
    sublist_file="${sample_prefix}_$((i+1)).list"  # Name of sublist file
    # Extract sublist
    sed -n "${start},${end}p" "$input_file" > "$PWD/input/pp5020/${formatted_date}/$sublist_file"  
done

echo $num_sublists
