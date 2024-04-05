#!/bin/bash

# Set path to CMSSW
source $HOME/setup_cmsenv.sh

# Initial parameters to run
EXEC_PATH=${HOME}/soft/jetAnalysis/processing
cd $EXEC_PATH

files_per_job=2
#input_file_list=$HOME/filelists/pythiaHydjet2018_miniAODforest.txt
input_file_list=$HOME/filelists/test_list.txt

echo -e "Splitting input file list: ${input_file_list}"
n_sublists=$(./split_input_2_sublists.sh ${input_file_list} ${files_per_job})
echo -e "Input file list is splitted into ${n_sublists}"

for ((jobId = 1; jobId <= $n_sublists; jobId++)); do
cat <<EOF >condor/sub/job_jetAna_$jobId.sub
universe = vanilla

executable = run_jetAna.sh
+JobFlavour           = "longlunch"
getenv     = True
requirements =((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"

arguments             = input/inputfile_$jobId.list output/jetAna/oJetAna_$jobId.root

output                = condor/log/jetAna.${jobId}.out
error                 = condor/log/jetAna.${jobId}.err
log                   = condor/log/jetAna.${jobId}.log

queue 

EOF

    condor_submit condor/sub/job_jetAna_$jobId.sub
done

