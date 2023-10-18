// C++ headers
#include <iostream>

// Jet analysis headers
#include "Manager.h"
#include "ForestAODReader.h"
#include "JetAnalysis.h"
#include "BasicHistoManager.h"
#include "EventCut.h"
#include "JetCut.h"
#include "TFile.h"

//________________
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {
    /* code */

    // Set default values for arguments
    TString inFileName = "../../../data/HiForestAOD_PbPbMC2018skim_10.root";
    TString oFileName = "oTestSimpleReadForest.root";
    Long64_t nEventsToRead = 500;

    // Read input argument list 
    if (argc > 1) inFileName = argv[1];
    if (argc > 2) oFileName = argv[2];

    Manager *manager = new Manager();
    EventCut *eventCut = new EventCut();
    eventCut->setVz(-20., 20.);
    JetCut *jetCut = new JetCut();
    jetCut->setMustHaveGenMathing();
    jetCut->setRecoPt(50., 1500.);
    //jetCut->setVerbose();
    ForestAODReader *reader = new ForestAODReader(inFileName);
    reader->setIsMc(kTRUE);
    reader->usePartFlowJetBranch();
    //reader->useCaloJetBranch();
    reader->setCollidingSystem("PbPb");
    reader->setCollidingEnergy(5020);
    reader->setYearOfDataTaking(2018);
    reader->setEventCut(eventCut);
    reader->setJetCut(jetCut);
    reader->fixJetArrays();
    manager->setEventReader(reader);

    JetAnalysis *analysis = new JetAnalysis();
    BasicHistoManager *hm = new BasicHistoManager();
    hm->setIsMc(kTRUE);
    hm->init(kTRUE); // kTRUE stands up for use MC; need to FIX
    analysis->addHistoManager(hm);
    manager->addAnalysis(analysis);

    manager->init();
    manager->performAnalysis();
    manager->finish();

    TFile* oFile = new TFile(oFileName, "recreate");
    hm->writeOutput();
    oFile->Close();
    
    return 0;
}
