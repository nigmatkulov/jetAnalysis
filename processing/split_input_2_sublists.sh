#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <N>"
    exit 1
fi

input_file="$1"  # Input file name
N="$2"           # Number of entries for each sublist

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

# Create sublists
for ((i = 0; i < num_sublists; i++)); do
    start=$((i * N + 1))  # Calculate start line number for current sublist
    end=$((start + N - 1))  # Calculate end line number for current sublist
    sublist_file="inputfile_$((i+1)).list"  # Name of sublist file
    sed -n "${start},${end}p" "$input_file" > "$PWD/input/$sublist_file"  # Extract sublist
done

echo $num_sublists

