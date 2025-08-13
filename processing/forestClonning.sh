#!/bin/bash

# Source the CMSSW environment
source $HOME/setup_cmsenv.sh

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input ROOT filename> <output directory>"
    exit 1
fi

# Get the input arguments
inputFileName=$1
outputDirectory=$2

# Check if the input file exists
if [ ! -f "$inputFileName" ]; then
    echo "Error: Input file '$inputFileName' not found!"
    exit 1
fi

# Check if the output directory exists
if [ ! -d "$outputDirectory" ]; then
    echo "Error: Output directory '$outputDirectory' not found! Creating the one..."
    mkdir -p $outputDirectory
    echo "Output directory '$outputDirectory' created!"  
fi

# Call the cloneFile.C macro using ROOT
if [ $? -eq 0 ]; then
    root -l -b -q "~/soft/jetAnalysis/macro/cloneForest.C(\"$inputFileName\", \"$outputDirectory\")"
else
    echo "Error: Failed to execute macro"
    exit 1
fi
