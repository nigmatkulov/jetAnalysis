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
jeuSyst=$7
jerSyst=$8

echo -e "Input file list             : ${input_file_list}"
echo -e "Output file name            : ${output_file_name}"
echo -e "Is MC                       : ${is_mc}"
echo -e "isPbGoingDir (RunB or RunD) : ${is_Pbgoing}"
echo -e "ptHatLow                    : ${pt_hat_low}"
echo -e "ptHatHi                     : ${pt_hat_hi}"
echo -e "JEU syst         : ${jeuSyst}"
echo -e "JER syst         : ${jerSyst}"

# Run jetAna
if [ "$is_mc" -eq 1 ]; then
    if [ "$is_Pbgoing" -eq 1 ]; then
        ../build/dijetAna_pPb5020 ${input_file_list} /eos/user/g/gnigmatk/ana/pPb5020/embedding/RunD/${output_file_name} ${is_mc} ${is_Pbgoing} ${pt_hat_low} ${pt_hat_hi} ${jeuSyst} ${jerSyst}
		else
        ../build/dijetAna_pPb5020 ${input_file_list} /eos/user/g/gnigmatk/ana/pPb5020/embedding/RunB/${output_file_name} ${is_mc} ${is_Pbgoing} ${pt_hat_low} ${pt_hat_hi} ${jeuSyst} ${jerSyst}
    fi
else
    if [ "$is_Pbgoing" -eq 1 ]; then
        ../build/dijetAna_pPb5020 ${input_file_list} /eos/user/g/gnigmatk/ana/pPb5020/exp/RunD/${output_file_name} ${is_mc} 0 ${pt_hat_low} ${pt_hat_hi} ${jeuSyst} ${jerSyst}
    else
        ../build/dijetAna_pPb5020 ${input_file_list} /eos/user/g/gnigmatk/ana/pPb5020/exp/RunB/${output_file_name} ${is_mc} 0 ${pt_hat_low} ${pt_hat_hi} ${jeuSyst} ${jerSyst}
    fi
fi

echo -e "Data processing of ${input_file_list} is finished"
