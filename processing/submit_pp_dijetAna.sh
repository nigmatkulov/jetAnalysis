#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

# Date of submission
formatted_date=$(date +"%Y%m%d")

# Read input parameters
# First parameters tells the dataset name: DATA_MB, DATA_HM185, DATA_HM250, DATA_PAEGJet
#sample_name=PYTHIA
sample_name=DATA

# Beam direction
is_Pbgoing=1
if [ "$is_Pbgoing" -eq 1 ]; then
    direction=Pbgoing
else
    direction=pgoing
fi

# Generate path to the inputfile list
if [ "$sample_name" == "DATA" ]; then
    sample_prefix="ppData2017_jet60or80Triggers"
    input_file_list="${EXEC_PATH}/filelists/pp5020/DATA/${sample_prefix}.list"
    
elif [ "$sample_name" == "PYTHIA" ]; then
    sample_prefix="Pythia_pThat15_full"
    input_file_list="${EXEC_PATH}/filelists/pp5020/PYTHIA/${sample_prefix}.list"
fi

# Specify number of files per list to split
files_per_job=50
#input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_dijet_input.sh ${input_file_list} ${files_per_job} ${sample_name} ${direction} ${pd_number})
echo -e "Input file list is splitted into ${n_sublists}"

if [ ! -d "condor/sub/ppp5020/${formatted_date}" ]; then
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
    mkdir -p "condor/log/pPb8160/${formatted_date}"
    if [ $? -eq 0 ]; then
        echo "Directory created successfully."
    else
        echo "Failed to create directory. Terminating"
        exit 1
    fi
fi

cat <<EOF > condor/sub/pp5020/${formatted_date}/pp5020_${sample_prefix}.sub
universe = vanilla
executable = ${EXEC_PATH}/run_dijetAna.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
EOF


for ((jobId = 1; jobId <= $n_sublists; jobId++)); do
    cat <<EOF >>condor/sub/pp5020/${formatted_date}/pp5020_${sample_prefix}.sub
arguments             = input/pp5020/${formatted_date}/${sample_prefix}_$jobId.list ${sample_prefix}_pp5020_$jobId.root 0 ${is_Pbgoing} 0 15000
output                = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.out
error                 = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.err
log                   = condor/log/pp5020/${formatted_date}/${sample_prefix}_$jobId.log
queue 

EOF

done

condor_submit condor/sub/pPb8160/${formatted_date}/pPb8160_${sample_prefix}.sub


