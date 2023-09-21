// ROOT headers
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TDirectoryFile.h"
#include "TFile.h"

// C++ headers
#include <cmath>
#include <string>
#include <math.h>
#include <cstring>
#include <iostream>

// Jet analysis headers
#include "JetAnalysis.h"

ClassImp(JetAnalysis)

//________________
JetAnalysis::JetAnalysis(TString inFileName, 
                         TString oFileName, 
                         Long64_t nEventsToRead) :
                         TObject(),
                         fInFileName(inFileName),
                         fOutFileName(oFileName.Data()),
                         fEvent2Read(nEventsToRead),
                         fIsMc(kFALSE),
                         fHltTree(nullptr),
                         fEventTree(nullptr),
                         fJetTree(nullptr),
                         fTrkTree(nullptr) {

    std::cout << "Creating and instance of JetAnalysis";
    std::cout << "\t[DONE]\n";
}

//________________
JetAnalysis::~JetAnalysis() {
    /* Empty */
    Clear();
}

//________________
void JetAnalysis::init() {
    // Initialize analysis
    std::cout << "JetAnalysis::init" << std::endl;
}

//________________
void JetAnalysis::processData() {
    // Perform the analysis
    std::cout << "JetAnalysis::processData" << std::endl;
}

//________________
void JetAnalysis::finish() {
    // Save data and close files
    std::cout << "JetAnalysis::finfish" << std::endl;
}

//________________
void JetAnalysis::setupInput() {
    
}

//_________________
Bool_t JetAnalysis::fileExists() {
    
}