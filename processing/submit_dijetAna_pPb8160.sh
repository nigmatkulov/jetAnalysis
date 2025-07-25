#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Trigger case
triggerId=0 # 0 - MB, 1 - jet60, 2 - jet80, 3 - jet100

sample_name=DATA_MB
if [ $triggerId -ne 0 ]; then
    sample_name=DATA_PAEGJet
fi

# Beam direction
is_Pbgoing=1
if [ "$is_Pbgoing" -eq 1 ]; then
    direction=Pbgoing
else
    direction=pgoing
fi


# Dataset number
pd_number=$1

# JEU systematics: 0 - default, -1 - JEU-, 1 - JEU+
jeuSyst=0 
# JER systematics: 0 - default, -1 - JER-, 1 - JER+, other - no extra smearing is applied (pure JEC)
jerSyst=-99
# RecoJet selection method: 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId
recoJetSelMethod=1

# Generate path to the inputfile list
if [ "$sample_name" == "DATA_MB" ]; then
    sample_prefix="MB_PD${pd_number}_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb8160/DATA_MB/${direction}/${sample_prefix}.txt"
    
elif [ "$sample_name" == "DATA_HM185" ]; then
    sample_prefix="HM185_PD${pd_number}_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb8160/DATA_HM185/${direction}/${sample_prefix}.txt"
    
elif [ "$sample_name" == "DATA_HM250" ]; then
    sample_prefix="HM250_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb8160/DATA_HM250/${direction}/${sample_prefix}.txt"
else
    sample_prefix="PAEGJet_${direction}"
    input_file_list="${EXEC_PATH}/filelists/pPb8160/DATA_PAEGJet/${direction}/${sample_prefix}.txt"
fi

# Specify number of files per list to split
files_per_job=50
#input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_pPb8160_dataset.sh ${input_file_list} ${files_per_job} ${sample_name} ${direction} ${pd_number})
echo -e "Input file list is splitted into ${n_sublists}"

if [ ! -d "condor/sub/pPb8160/${formatted_date}" ]; then
    echo "Directory 'condor/sub/pPb8160/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/sub/pPb8160/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

if [ ! -d "condor/log/pPb8160/${formatted_date}" ]; then
    echo "Directory 'condor/log/pPb8160/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/log/pPb8160/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

cat <<EOF > condor/sub/pPb8160/${formatted_date}/pPb8160_${sample_prefix}.sub
universe = vanilla
executable = ${EXEC_PATH}/run_dijetAna_pPb8160.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
EOF


for ((jobId = 1; jobId <= $n_sublists; jobId++)); do
    cat <<EOF >>condor/sub/pPb8160/${formatted_date}/pPb8160_${sample_prefix}.sub
arguments             = input/pPb8160/${formatted_date}/${sample_prefix}_$jobId.list ${sample_prefix}_pPb8160_$jobId.root 0 ${is_Pbgoing} 0 15000 ${jeuSyst} ${jerSyst} ${triggerId} ${recoJetSelMethod}
output                = condor/log/pPb8160/${formatted_date}/${sample_prefix}_$jobId.out
error                 = condor/log/pPb8160/${formatted_date}/${sample_prefix}_$jobId.err
log                   = condor/log/pPb8160/${formatted_date}/${sample_prefix}_$jobId.log
queue 

EOF

done

condor_submit condor/sub/pPb8160/${formatted_date}/pPb8160_${sample_prefix}.sub


