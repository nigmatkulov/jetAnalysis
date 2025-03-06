// ROOT headers
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

// C++ headers
#include <iostream>
#include <vector>
#include <string>

// Function to clone a subset of trees from a ROOT file
void cloneForest(const char* inputFileName, const char *outputDirectory = "./") {
    // Open the original ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file!" << std::endl;
        return;
    }

    TString outputFileName = TString(outputDirectory);
    outputFileName.Append("/").Append(TString(inputFileName).ReplaceAll(".root", "_trunc.root"));
    // Create a new ROOT file
    TFile* outputFile = new TFile(outputFileName.Data(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file!" << std::endl;
        inputFile->Close();
        return;
    }

    // List of directories and their corresponding trees to copy
    std::vector<std::pair<std::string, std::string>> directoriesToCopy = {
        {"hiEvtAnalyzer", "HiTree"},
        {"hltanalysis", "HltTree"},
        {"skimanalysis", "HltTree"},
        {"ak4PFJetAnalyzer", "t"} /*,
        {"akCs4PFJetAnalyzer", "t"},
        {"ak3PFJetAnalyzer", "t"}
        */
    };

    // Loop through selected directories
    for (const auto& entry : directoriesToCopy) {
        const std::string& dirName = entry.first;
        const std::string& treeName = entry.second;

        // Get directory from input file
        TDirectory* inputDir = (TDirectory*)inputFile->Get(dirName.c_str());
        if (!inputDir) {
            std::cerr << "Warning: Directory " << dirName << " not found!" << std::endl;
            continue;
        }

        // Create the same directory in the output file
        outputFile->mkdir(dirName.c_str());
        outputFile->cd(dirName.c_str());  // Move into the new directory

        // Get the tree in the directory
        TTree* tree = (TTree*)inputDir->Get(treeName.c_str());
        if (!tree) {
            std::cerr << "Warning: Tree " << treeName << " not found in directory " << dirName << std::endl;
            continue;
        }

        // Clone the tree into the new directory
        TTree* newTree = tree->CloneTree(-1);
        newTree->Write(); // Write the tree to the new file

        std::cout << "Copied directory: " << dirName << " with tree: " << treeName << std::endl;
    }

    // Close files
    outputFile->Close();
    inputFile->Close();

    std::cout << "New ROOT file created: " << outputFileName.Data() << std::endl;
}