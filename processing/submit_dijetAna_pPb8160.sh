#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

echo "Number of arguments passed: $#"
echo "Arguments: $@"

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Trigger case
trigger_id=1 # 0 - MB, 1 - jet60, 2 - jet80, 3 - jet100
trigger_name=MB
if [ $trigger_id -eq 1 ]; then
    trigger_name=Jet60
elif [ $trigger_id -eq 2 ]; then
    trigger_name=Jet80
elif [ $trigger_id -eq 3 ]; then
    trigger_name=Jet100
fi

# Beam direction
is_pbgoing=1
if [ "$is_pbgoing" -eq 1 ]; then
    direction=Pbgoing
else
    direction=pgoing
fi

# Dataset number
if [ $# -eq 1 ]; then
    pd_number=$1 # PD number is defined for MB only. Pb-going 1-20, p-going 1-8
else
    pd_number=-1 
fi

# JEU systematics: 0 - default, -1 - JEU-, 1 - JEU+
jeu_syst=0 
# JER systematics: 0 - default, -1 - JER-, 1 - JER+, other - no extra smearing is applied (pure JEC)
jer_syst=-99
# RecoJet selection method: 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId
reco_jet_sel_method=2

# Specify number of files per list to split
files_per_job=50
#input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_pPb8160_dataset.sh -n ${files_per_job} -t ${trigger_id} -d ${direction} -p ${pd_number})
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

if [ "$trigger_id" -eq 0 ]; then
    prefix_name=${trigger_name}_PD${pd_number}_${direction}
else
    prefix_name=${trigger_name}_PAEG_${direction}
fi

# Loop over the sublists
for ((job_id = 1; job_id <= $n_sublists; job_id++)); do
    # Get filename with jobId attached
    prefix_name_with_job_id=${prefix_name}_${job_id}
    # Create the submission file
    cat <<EOF > condor/sub/pPb8160/${formatted_date}/pPb8160_${prefix_name_with_job_id}.sub
universe = vanilla
executable = ${EXEC_PATH}/run_dijetAna_pPb8160.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"

arguments             = input/pPb8160/${formatted_date}/${prefix_name_with_job_id}.list ${prefix_name_with_job_id}.root 0 ${is_pbgoing} 0 15000 ${jeu_syst} ${jer_syst} ${trigger_id} ${reco_jet_sel_method}
output                = condor/log/pPb8160/${formatted_date}/${prefix_name_with_job_id}.out
error                 = condor/log/pPb8160/${formatted_date}/${prefix_name_with_job_id}.err
log                   = condor/log/pPb8160/${formatted_date}/${prefix_name_with_job_id}.log
queue

EOF
    
    # Submit the job
    condor_submit condor/sub/pPb8160/${formatted_date}/pPb8160_${prefix_name_with_job_id}.sub
done

echo "All jobs for ${prefix_name} submitted."
