#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Read input parameters
# First parameters tells the dataset name: 0 - DATA_MB, 1 - DATA_PAEGJet
#sample_name=DATA_MB
sample_name=1

# JEU systematics: 0 - default, -1 - JEU-, 1 - JEU+
jeuSyst=0 
# JER systematics: 0 - default, -1 - JER-, 1 - JER+
jerSyst=0

# Beam direction: 0 - RunB, 1 - RunD
is_Pbgoing=1
if [ "$is_Pbgoing" -eq 1 ]; then
    direction=RunD
else
    direction=RunB
fi

# Generate path to the inputfile list
if [ "$sample_name" -eq 0 ]; then
    sample_prefix="MB_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb5020/Data_MB/${direction}/MB_${direction}5TeV_allPDs.txt"
else
    sample_prefix="PAEGJet_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb5020/Data_PAEGJet/${direction}/PAEGJet_${direction}5TeV.txt"
fi

# Specify number of files per list to split
files_per_job=100
#input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_pPb5020_dataset.sh ${input_file_list} ${files_per_job} ${sample_name} ${is_Pbgoing})
echo -e "Input file list is splitted into ${n_sublists}"

if [ ! -d "condor/sub/pPb5020/${formatted_date}" ]; then
    echo "Directory 'condor/sub/pPb5020/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/sub/pPb5020/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

if [ ! -d "condor/log/pPb5020/${formatted_date}" ]; then
    echo "Directory 'condor/log/pPb5020/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/log/pPb5020/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

cat <<EOF > condor/sub/pPb5020/${formatted_date}/pPb5020_${sample_prefix}.sub
universe = vanilla
executable = ${EXEC_PATH}/run_dijetAna_pPb5020.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
EOF


for ((jobId = 1; jobId <= $n_sublists; jobId++)); do
    cat <<EOF >>condor/sub/pPb5020/${formatted_date}/pPb5020_${sample_prefix}.sub
arguments             = input/pPb5020/${formatted_date}/${sample_prefix}_$jobId.list ${sample_prefix}_pPb5020_$jobId.root 0 ${is_Pbgoing} 0 15000 ${jeuSyst} ${jerSyst}
output                = condor/log/pPb5020/${formatted_date}/${sample_prefix}_$jobId.out
error                 = condor/log/pPb5020/${formatted_date}/${sample_prefix}_$jobId.err
log                   = condor/log/pPb5020/${formatted_date}/${sample_prefix}_$jobId.log
queue 

EOF

done

condor_submit condor/sub/pPb5020/${formatted_date}/pPb5020_${sample_prefix}.sub


