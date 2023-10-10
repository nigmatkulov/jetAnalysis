// ROOT headers
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"

// C++ headers
#include <iostream>

// Jet analysis headers
#include "JetAnalysis.h"

ClassImp(JetAnalysis)

//________________
JetAnalysis::JetAnalysis() : BaseAnalysis(), fDebug(kFALSE) {
    /* Empty */
}

//________________
JetAnalysis::~JetAnalysis() {
    /* Empty */
}

//________________
void JetAnalysis::init() {
    // Initialize analysis
    std::cout << "JetAnalysis::init" << std::endl;
}

//________________
void JetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    std::cout << "JetAnalysis::processEvent" << std::endl;
    //std::cout << "vz: " << event->vz() << std::endl;
}

//________________
void JetAnalysis::finish() {
    // Save data and close files
    std::cout << "JetAnalysis::finish" << std::endl;
}

//________________
void JetAnalysis::report() {
    // Force to report everyone
}

//________________
TList* JetAnalysis::getOutputList() {
    TList *outputList = new TList();

    // Add list of settings for cuts

    return outputList;
}