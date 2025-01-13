#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Read input parameters
# First parameters tells the dataset name: 0 - DATA, 1 - PYTHIA
sample_name=0

# JEU systematics: 0 - default, -1 - JEU-, 1 - JEU+
jeuSyst=0 
# JER systematics: 0 - default, -1 - JER-, 1 - JER+
jerSyst=0

# Generate path to the inputfile list
if [ "$sample_name" -eq 0 ]; then
    sample_prefix="data"
    input_file_list="${EXEC_PATH}/filelists/pp5020/DATA/ppData2017_jet60or80Triggers.list"
else
    sample_prefix="pythia"
    input_file_list="${EXEC_PATH}/filelists/pp5020/PYTHIA/Pythia_pThat15.list"
fi

# Specify number of files per list to split
files_per_job=100
#input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_pp5020_dataset.sh ${input_file_list} ${files_per_job} ${sample_name})
echo -e "Input file list is splitted into ${n_sublists}"

if [ ! -d "condor/sub/pp5020/${formatted_date}" ]; then
    echo "Directory 'condor/sub/pp5020/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/sub/pp5020/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

if [ ! -d "condor/log/pp5020/${formatted_date}" ]; then
    echo "Directory 'condor/log/pp5020/${formatted_date}' does not exist. Creating..."
    mkdir -p "condor/log/pp5020/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

cat <<EOF > condor/sub/pp5020/${formatted_date}/${sample_prefix}.sub
universe = vanilla
executable = ${EXEC_PATH}/run_dijetAna_pp5020.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
EOF


for ((jobId = 1; jobId <= $n_sublists; jobId++)); do
    cat <<EOF >>condor/sub/pp5020/${formatted_date}/${sample_prefix}.sub
arguments             = input/pp5020/${formatted_date}/${sample_prefix}_$jobId.list ${sample_prefix}_$jobId.root ${sample_name} 0 -100000000 100000000 ${jeuSyst} ${jerSyst}
output                = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.out
error                 = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.err
log                   = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.log
queue 

EOF

done

condor_submit condor/sub/pp5020/${formatted_date}/${sample_prefix}.sub


