#!/bin/bash

# Source the CMSSW environment
source $HOME/setup_cmsenv.sh

EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

usage() {
    echo "Usage: $0 -t <trigger> -d <direction> [-p <PD>] -o <output_dir>"
    echo "  -t <trigger>: 0 (MB) or 1 (PAEGJet)"
    echo "  -d <direction>: 1 (Pbgoing) or 0 (pgoing)"
    echo "  -p <PD>: Physics Dataset name (required for MB trigger)"
    echo "  -o <output_dir>: Path to output directory"
    exit 1
}

while getopts ":t:d:p:o:" opt; do
    case $opt in
        t) trigger="$OPTARG" ;;
        d) direction="$OPTARG" ;;
        p) PD="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$trigger" || -z "$direction" || -z "$output_dir" ]]; then
    usage
fi

if [[ "$trigger" != "0" && "$trigger" != "1" ]]; then
    echo "Error: trigger must be 0 (MB) or 1 (PAEGJet)"
    usage
fi

if [[ "$direction" != "0" && "$direction" != "1" ]]; then
    echo "Error: direction must be 0 (pgoing) or 1 (Pbgoing)"
    usage
fi

if [[ "$trigger" == "0" && -z "$PD" ]]; then
    echo "Error: PD must be specified for MB trigger"
    usage
fi

if [[ "$direction" == "0" ]]; then
    direction_name="pgoing"
else
    direction_name="Pbgoing"
fi

if [[ "$trigger" == "0" ]]; then
    trigger_name="MB${PD}/$direction_name"
else
    trigger_name="PAEGJet/$direction_name"
fi


# Set the input file based on the trigger and direction
if [[ "$trigger" == "0" ]]; then
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_MB/${direction_name}/MB_PD${PD}_${direction_name}.txt"
else
    input_file="${EXEC_PATH}/filelists/pPb8160/DATA_PAEGJet/${direction_name}/PAEGJet_${direction_name}.txt"
fi

echo "Using input file: $input_file"

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Check if output directory exists, create if not
if [[ ! -d "$output_dir" ]]; then
    mkdir -p "$output_dir"
    echo "Created directory: $output_dir"
fi

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Check if directory in the condor/sub/$formatted_date/foresting directory exists
path_2_sub_files="$EXEC_PATH/condor/sub/pPb8160/$formatted_date/foresting/$trigger_name"
if [[ ! -d "$path_2_sub_files" ]]; then
    mkdir -p "$path_2_sub_files"
    echo "Created directory: $path_2_sub_files"
fi

# Check if directory to log, err, and out files exists
path_2_log_files="condor/log/pPb8160/$formatted_date/foresting/$trigger_name"
if [[ ! -d "$EXEC_PATH/$path_2_log_files" ]]; then
    mkdir -p "$EXEC_PATH/$path_2_log_files"
    echo "Created directory: $EXEC_PATH/$path_2_log_files"
fi

echo "Start input file processing"
n_files=0
# Parse input file line by line
while IFS= read -r line; do
    n_files=$((n_files + 1))
    # Process each line
    if [[ -z "$line" ]]; then
        continue  # Skip empty lines
    fi
    input_file_basename=$(basename "$line")
    input_file_basename="${input_file_basename%.*}"
    sub_file="${path_2_sub_files}/${input_file_basename}.sub"
    cat <<EOF > $sub_file
universe = vanilla
executable = ${EXEC_PATH}/forestClonning.sh
+JobFlavour           = "microcentury"
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
getenv     = True
RequestCpus = 1
transfer_input_files  = ${EXEC_PATH}/voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"

arguments = $line $output_dir
output = ${path_2_log_files}/${input_file_basename}.out
log = ${path_2_log_files}/${input_file_basename}.log
error = ${path_2_log_files}/${input_file_basename}.err
queue
EOF
    
    # condor_submit $sub_file
done < $input_file

echo "Submitted $n_files jobs for forest cloning."
