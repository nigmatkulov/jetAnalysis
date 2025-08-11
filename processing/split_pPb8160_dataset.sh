#!/bin/bash

EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

n_files_per_sublist=50 # Number of files per sublist
trigger_id=0           # 0 - MB, 1 - jet60, 2 - jet80, 3 - jet100
direction="Pbgoing"    # Pbgoing or pgoing
pd_number=-1            # PD number is defined for MB only. Pb-going 1-20, p-going 1-8          

while getopts "n:t:d:p:" opt; do
    case $opt in
        n) n_files_per_sublist="$OPTARG" ;;  # -n <n_files_per_sublist>
        t) trigger_id="$OPTARG" ;;           # -t <trigger_id>
        d) direction="$OPTARG" ;;            # -d <direction>
        p) pd_number="$OPTARG" ;;            # -p <pd_number>
        # Add more options below as needed, e.g.:
        # x) new_var="$OPTARG" ;;             # -x <new_var>
        \?) echo "Usage: $0 [-n n_files_per_sublist] [-t trigger_id] [-d direction] [-p pd_number] [not supported other options]" >&2
                exit 1 ;;
    esac
done

echo "Arguments passed to the script: $@"
echo "n_files_per_sublist: $n_files_per_sublist"
echo "trigger_id: $trigger_id"
echo "direction: $direction"
echo "pd_number: $pd_number"

if [[ "$trigger_id" -eq 0 && "$pd_number" -eq -1 ]]; then
    echo "Error: PD number must be specified for MB trigger."
    exit 1
fi

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Trigger name
trigger_name="MB"
if [ "$trigger_id" -eq 1 ]; then
    trigger_name="Jet60"
elif [ "$trigger_id" -eq 2 ]; then
    trigger_name="Jet80"
elif [ "$trigger_id" -eq 3 ]; then
    trigger_name="Jet100"
fi

# Sample name
sample_name="DATA_MB"
if [ "$trigger_id" -ne 0 ]; then
    sample_name="DATA_PAEGJet"
fi

# Generate path to the inputfile list (IMPORTANT:HM triggers omitted)
if [ "$sample_name" == "DATA_MB" ]; then
    sample_prefix="MB_PD${pd_number}_${direction}"
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_MB/${direction}/${sample_prefix}.txt"
elif [ "$sample_name" == "DATA_HM185" ]; then
    sample_prefix="HM185_PD${pd_number}_${direction}"
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_HM185/${direction}/${sample_prefix}.txt"
    
elif [ "$sample_name" == "DATA_HM250" ]; then
    sample_prefix="HM250_${direction}"
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_HM250/${direction}/${sample_prefix}.txt"
else
    sample_prefix="PAEGJet_${direction}"
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_PAEGJet/${direction}/${sample_prefix}.txt"
fi


# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Input file '$input_file' not found."
    exit 1
fi

# Count the total number of lines in the input file
total_lines=$(wc -l < "$input_file")

# Calculate the number of sublists needed
num_sublists=$((total_lines / ${n_files_per_sublist}))
if [ $((total_lines % ${n_files_per_sublist})) -ne 0 ]; then
    ((num_sublists++))
fi

if [ ! -d "$PWD/input/pPb8160/${formatted_date}" ]; then
    #echo "Directory '$PWD/input/pPb8160/${formatted_date}' does not exist. Creating..."
    mkdir -p "$PWD/input/pPb8160/${formatted_date}"
fi

# Create sublists
for ((i = 0; i < num_sublists; i++)); do
    start=$((i * n_files_per_sublist + 1))  # Calculate start line number for current sublist
    end=$((start + n_files_per_sublist - 1))  # Calculate end line number for current sublist
    if [ "$trigger_id" -eq 0 ]; then
        sublist_file="${trigger_name}_PD${pd_number}_${direction}_$((i+1)).list"  # Name of sublist file for PAEGJet
    else
        sublist_file="${trigger_name}_${direction}_$((i+1)).list"  # Name of sublist file
    fi
    # Extract sublist
    sed -n "${start},${end}p" "$input_file" > "$PWD/input/pPb8160/${formatted_date}/$sublist_file"  
done

echo $num_sublists
