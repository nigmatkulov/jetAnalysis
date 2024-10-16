#!/bin/bash

# Load CMS environment
source ${HOME}/setup_cmsenv.sh

# Add path to the PDF sets

# Where the executable and the library are stored
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

input_file_list=$1
output_file_name=$2
is_mc=$3
is_Pbgoing=$4
pt_hat_low=$5
pt_hat_hi=$6

echo -e "Input file list  : ${input_file_list}"
echo -e "Output file name : ${output_file_name}"
echo -e "Is MC            : ${is_mc}"
echo -e "isPbGoingDir     : ${is_Pbgoing}"
echo -e "ptHatLow         : ${pt_hat_low}"
echo -e "ptHatHi          : ${pt_hat_hi}"

# Run jetAna
../build/dijetAna ${input_file_list} ${output_file_name} ${is_mc} ${is_Pbgoing} ${pt_hat_low} ${pt_hat_hi}


echo -e "Data processing of thes ${input_file_list} is finished"
