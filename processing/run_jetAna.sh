#!/bin/bash

# Load CMS environment
source ${HOME}/setup_cmsenv.sh

# Add path to the PDF sets

# Where the executable and the library are stored
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

input_file_list=$1
output_file_name=$2

echo -e "Input file list: ${input_file_list}"
echo -e "Output file name: ${output_file_name}"

# Run jetAna
../build/jetAna ${input_file_list} ${output_file_name}

echo -e "Data processing of thes ${input_file_list} is finished"
